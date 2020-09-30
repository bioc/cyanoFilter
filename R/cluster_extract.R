#' extract clusters based on supplied cluster indicator
#'
#' @param flowfile flowframe containing cluster indicators as well
#' @param cluster_var column name in expression matrix containing the cluter indicators, cannot be NULL.
#' @param cluster_val cluster number, cannot be NULL.
#'
#'
#'@examples
#'
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'                              mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = cells_nonmargin$reducedflowframe,
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' fin <- phyto_filter(flowfile = cells_nodebris$syn,
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' cluster_extract(flowfile = fin$flowframe_proportion,
#'     cluster_var = "Clusters",
#'     cluster_val = 1)
#'
#' @export cluster_extract



cluster_extract <- function(flowfile,
                            cluster_var = "Clusters",
                            cluster_val = NULL) {


  if(is.null(cluster_var) | is.null(cluster_val)) {

    stop("Either of cluster_var or cluster_val is empty")

  } else {

          return(flowfile[flowfile@exprs[, cluster_var] %in% cluster_val, ])

  }


}



#' takes a flowframe, name of cluster column and extracts part of flowframe that makes up proportion.
#'
#' @param flowfile flowframe after debris are removed.
#' @param cluster_var column name in expression matrix containing the cluter indicators
#' @param proportion value between 0 and 1 indicating percentage of the total particles wanted
#'
#' @export cluster_extractp


cluster_extractp <- function(flowfile, cluster_var = "Clusters", proportion = 0.80) {

  if(is.null(proportion)) {

    proportion <- 1
    message("Proportion set to 1, hence attempt to recover all particles measured")

  }

  if(is.null(cluster_var)) {

    stop("Error: Cluster_var is not supplied")

  } else {

    cluster_n <- data.frame(table(flowfile@exprs[,cluster_var]))
    names(cluster_n) <- c("Cluster", "Frequency")
    cluster_n$Percentage <- cluster_n$Frequency/sum(cluster_n$Frequency)
    cluster_n <- cluster_n[order(cluster_n$Percentage, decreasing = TRUE),]
    cluster_n$CPercentage <- cumsum(cluster_n$Percentage)
    c1 <- cluster_n$Cluster[which(cluster_n$CPercentage < proportion)]
    c2 <- cluster_n$Cluster[which(cluster_n$CPercentage >= proportion)[1]]
    needed_clusters <- c(c1, c2)
    needed_flowfile <- flowfile[flowfile@exprs[, cluster_var] %in% needed_clusters, ]

    return(list(particles_per_cluster = cluster_n,
                clusters_proportion = needed_clusters,
                flowfile_proportion = needed_flowfile)
           )

  }


}
