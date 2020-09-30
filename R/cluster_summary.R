#' takes a flowframes, a vector of channels, cluster indicator and return desired summaries per cluster
#'
#' @param object An object of class cyanoFilter to be summarised.
#' @param channels channels whose summaries are to be computed
#' @param cluster_var column name in expression matrix containing the cluter indicators
#' @param summary summary statistic of interest. Only mean and variance-covariance matrix supported at the moment.
#' @param ... other arguments. Not used at the moment
#'
#' @examples
#'
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = cells_nonmargin$reducedflowframe,
#'                            ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' fin <- phyto_filter(flowfile = cells_nodebris$syn,
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' summary(object = fin$flowframe_proportion,
#'         channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'         cluster_var = "Clusters",
#'         summary = 'mean')
#'
#' @importFrom stats aggregate median sd
#' @export
 summary <- function(object, channels = NULL,
                     cluster_var = "Clusters",
                     summary = c("mean", "median" , "cov", "n"), ...) {

   UseMethod('summary')

}

#' @export
summary.cyanoFilter <- function(object, channels = NULL,
                            cluster_var = "Clusters",
                            summary = c("mean", "median" , "cov", "n"), ...) {

  if(sum(summary %in% c("mean", "median", "cov", "n")) == 0) {

    stop("wrong summary option supplied")

  } else if(class(object)[2] == 'phytoFilter' & (is.null(channels) | is.null(cluster_var)) ) {

    stop("You must supply channels and cluster_var for objects of class phytoFilter")
  }

  if(class(object)[2] == 'phytoFilter') {

    sm <- match.arg(summary)

    ff <- object$flowframe_proportion

    if(sm == "mean") {

      #calculate mean per cluster
      mns <- aggregate(ff@exprs[, channels],
                       list(ff@exprs[, cluster_var]),
                       mean, na.rm = TRUE)
      names(mns) <- c(cluster_var, channels)
      return(mns)

    } else if(sm == "median") {

      #calculate mean per cluster
      mns <- aggregate(ff@exprs[, channels],
                       list(ff@exprs[, cluster_var]),
                       median, na.rm = TRUE)
      names(mns) <- c(cluster_var, channels)
      return(mns)

    } else if(sm == "n") {

      #calculate mean per cluster
      mns <- aggregate(ff@exprs[, channels],
                       list(ff@exprs[, cluster_var]),
                       length)[, 1:2]
      names(mns) <- c(cluster_var, "Number_of_particles")
      return(mns)

    } else if(sm == "cov") {

      #calculate variance per cluster
      uqs <- unique(ff@exprs[, cluster_var])
      vars_list <- lapply(uqs, function(i) {

        cov(ff@exprs[ff@exprs[, cluster_var] == i, channels])

      })
      names(vars_list) <- paste("Clusters", uqs, sep = "_")
      return(vars_list)

    } else {

      #calculate mean per cluster
      mns <- aggregate(ff@exprs[, channels],
                       list(ff@exprs[, cluster_var]),
                       mean, na.rm = TRUE)
      names(mns) <- c("Clusters", channels)

      #calculate median per cluster
      mns2 <- aggregate(ff@exprs[, channels],
                       list(ff@exprs[, cluster_var]),
                       median, na.rm = TRUE)
      names(mns2) <- c("Clusters", channels)

      #calculate variance per cluster
      uqs <- unique(ff@exprs[, cluster_var])
      vars_list <- lapply(uqs, function(i) {

        cov(ff@exprs[ff@exprs[, cluster_var] == i, channels])

      })
      names(vars_list) <- paste("Clusters", uqs, sep = "_")
      return(list(means = mns, medians = mns2, variances = vars_list))
    }

  } else if(class(object)[2] == 'marginEvents') {

    ff <- object$reducedflowframe

    #calculate overall mean
    mns <- apply(ff@exprs[, channels], 2,
                     mean, na.rm = TRUE)
    #calculate overall mean
    mns2 <- apply(ff@exprs[, channels], 2,
                 median, na.rm = TRUE)
    #calculate overall variance covariance
    vars_list <- cov(ff@exprs[, channels])

    #Number of margin events
    nMargin <- object$N_margin
    #Number of nonNamrgin events
    nonMargin <- object$N_nonmargin

    return(list(means = mns, medians = mns2, variances = vars_list,
                marginCount = nMargin,
                nonMarginCount = nonMargin))


  } else if(class(object)[2] == 'Debris') {

    ff <- object$syn

    #calculate overall mean
    mns <- apply(ff@exprs[, channels], 2,
                 mean, na.rm = TRUE)
    #calculate overall mean
    mns2 <- apply(ff@exprs[, channels], 2,
                  median, na.rm = TRUE)
    #calculate overall variance covariance
    vars_list <- cov(ff@exprs[, channels])

    return(list(means = mns, medians = mns2, variances = vars_list))


  } else stop('object of wrong class supplied')




}


