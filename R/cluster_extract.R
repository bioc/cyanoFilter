#' extract clusters based on supplied cluster indicator
#'
#' @param flowfile flowframe containing cluster indicators as well
#' @param cluster_var column name in expression matrix containing the cluter 
#'                    indicators, cannot be NULL.
#' @param cluster_val cluster number, cannot be NULL.
#'
#' @return flowFrame containing the clusters
#' 
#'@examples
#'
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                              package = "cyanoFilter",
#'                              mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                                          c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellMargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' fin <- phytoFilter(flowfile = reducedFlowframe(cells_nonmargin),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' clusterExtract(flowfile = reducedFlowframe(fin),
#'     cluster_var = "Clusters",
#'     cluster_val = 1)
#'
#' @export clusterExtract



clusterExtract <- function(flowfile,
                            cluster_var = "Clusters",
                            cluster_val = NULL) {


  if(is.null(cluster_var) | is.null(cluster_val)) {

    stop("Either of cluster_var or cluster_val is empty")

  } else {

    return(flowfile[flowfile@exprs[, cluster_var] %in% cluster_val, ])

  }


}



#' takes a flowframe, name of cluster column and extracts part of flowframe 
#' that makes up proportion.
#'
#' @param flowfile flowframe after debris are removed.
#' @param cluster_var column name in expression matrix containing the 
#' cluter indicators
#' @param proportion value between 0 and 1 indicating percentage of the total 
#' particles wanted
#'
#' @return a list containing \itemize{
#' \item \strong{particles_per_cluster}
#' \item \strong{clusters_proportion}
#' \item \strong{flowfile_proportion}
#' }
#'@examples
#'
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                              package = "cyanoFilter",
#'                              mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                                          c('SSC.W', 'TIME'))
#' cells_nonmargin <- cyanoFilter::cellMargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' fin <- phytoFilter(flowfile = reducedFlowframe(cells_nonmargin),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' clusterExtractp(flowfile = reducedFlowframe(fin),
#'     cluster_var = "Clusters",
#'     proportion = 0.80)
#' @export clusterExtractp


clusterExtractp <- function(flowfile, cluster_var = "Clusters", 
  proportion = 1) {

  if(is.null(proportion)) {

    message("Proportion set to 1, hence attempt to recover all 
            particles measured")

  }

  if(is.null(cluster_var)) {

    stop("cluster_var is not supplied")

  } else {

    cluster_n <- data.frame(table(exprs(flowfile)[,cluster_var]))
    names(cluster_n) <- c("Cluster", "Frequency")
    cluster_n$Percentage <- cluster_n$Frequency/sum(cluster_n$Frequency)
    cluster_n <- cluster_n[order(cluster_n$Percentage, decreasing = TRUE),]
    cluster_n$CPercentage <- cumsum(cluster_n$Percentage)
    c1 <- cluster_n$Cluster[which(cluster_n$CPercentage < proportion)]
    c2 <- cluster_n$Cluster[which(cluster_n$CPercentage >= proportion)[1]]
    needed_clusters <- c(c1, c2)
    needed_flowfile <- flowfile[exprs(flowfile)[, cluster_var] %in% 
                                  needed_clusters, ]

    return(list(particles_per_cluster = cluster_n,
                clusters_proportion = needed_clusters,
                flowfile_proportion = needed_flowfile)
           )

  }


}
