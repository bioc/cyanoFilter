#' gates out and assign indicators to phytoplankton cells based on 
#' the expression of  measured cell complexity channels.
#'
#' @param flowfile flowframe after debris are removed.
#'
#' @param pig_channels flowcytometer channels measuring cell pigments.
#'
#' @param com_channels flowcytometer channels measuring cell complexity.
#'
#' @param ph maximum peak height to be ignored. This allows ignoring of tiny 
#'           peaks that could affect the gating process.
#'
#' @param proportion proportion of cell count to be returned.
#'
#'
#' @return object of class \link{PhytopFilter} containing; \itemize{
#' \item \strong{fullflowframe -} flowframe containing all phytoplankton cells 
#'                                with added columns indicating cluster
#' \item \strong{flowframe_proportion -} a part of fullflowframe containing 
#'                                       proportion of cell count.
#' \item \strong{clusters_proportion -} proportion of cells in each cluster
#' \item \strong{particles_per_cluster -} number of particles per cluster
#' \item \strong{Cluster_ind -} indicator for each cluster
#' \item \strong{gated_channels -} channels with multiple peaks
#' }
#'
#' @description This function takes in a flowframe with debris removed and 
#'              identifies the different phytoplankton cell population 
#'              based on cell pigmentation and/or
#'              complexity.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the 
#'          \emph{flowDensity} package to
#'          identify peaks and identify cut-off points between these peaks.
#'
#'@examples
#'  flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellMargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debrisNc(flowframe = reducedFlowframe(cells_nonmargin),
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' phytoFilter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#' 
#'
#'
#' @export phytoFilter


phytoFilter <- function(flowfile, pig_channels = NULL, 
                         com_channels = NULL,
                         ph = 0.05, proportion = 0.80) {

  if(is.null(pig_channels) & is.null(com_channels)) {

    stop("both pigment and cell complecity channels cannot be NULL. 
         Check your input")

  } else {

    
    #gate based on the pigments
    channels <- c(pig_channels, com_channels)
  }
  
  

  pigment_gating <- pigmentGate(flowfile = flowfile, pig_channels = channels,
                                 ph = ph)
  gen_togate <- pigment_gating$gated_channels
  fgate <- pigment_gating


  if(sum(is.na(gen_togate)) == length(gen_togate)) {
        # no channel was gated here

        full_flowframe <- newFlowframe(flowfile,
                                       group = rep(1, nrow(flowfile)), 
                                       togate = NULL)
        #plott1 <- ggpairsDens(full_flowframe, channels = channels, 
        #group = "Clusters")
        message("Only one cluster found across the supplied channels.")
        
        ret_result <- PhytopFilter(fullflowframe = full_flowframe,
                                  flowframe_proportion = full_flowframe,
                                  clusters_proportion = 0,
                                  particles_per_cluster = data.frame(),
                                  Cluster_ind = 0,
                                  gated_channels = "",
                                  channels = c(pig_channels, com_channels)
                                )
        return(ret_result)

    } else {

      clusters <- flowCore::exprs(fgate$full_flowframe)[, 
                                             paste(gen_togate[
                                               !is.na(gen_togate)], 
                                               "Cluster", sep = "_")]
      if(is.null(dim(clusters))) {

        is <- clusters

      } else {

        cluster_unique <- unique(clusters)

        is <- apply(clusters, 1, function(x) {

          ii <- c()
          for(i in seq_len(nrow(cluster_unique))) {

            if(sum(x == cluster_unique[i, ]) == ncol(cluster_unique)) 
              ii <- c(ii, i) else next
          }

          return(ii)
        })

      }

      ### formulate a master full flow file
      full_flowframe <- newFlowframe(flowfile, group = is, 
                                      togate = NULL)
      needed_proportion <- clusterExtractp(full_flowframe, 
                                            cluster_var = "Clusters",
                                            proportion = proportion)
      flowframe_proportion <- needed_proportion$flowfile_proportion
      clusters_proportion <- as.numeric(needed_proportion$clusters_proportion)
      particles_per_cluster <- needed_proportion$particles_per_cluster

      #plott1 <- ggpairsDens(flowframe_proportion, channels = channels, 
      #group = "Clusters")

      ret_result <- PhytopFilter(
        fullflowframe = full_flowframe,
        flowframe_proportion = flowframe_proportion,
        clusters_proportion = clusters_proportion,
        particles_per_cluster = particles_per_cluster,
        Cluster_ind = is,
        gated_channels = gen_togate[!is.na(gen_togate)],
        channels = c(pig_channels, com_channels)
      )

      #attr(ret_result, 'class') <- c('cyanoFilter', 'phytoFilter')
      #cyanoFilterPhyto(ret_result)
      #return(ret_result)
    }
}
