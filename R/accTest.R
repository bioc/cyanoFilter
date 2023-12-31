
#' tests the accuracy of several automated gating functions on monoculture
#' flow cytometry experiments.
#'
#' @param flowfile flowSet with each flowFrame being a phytoplankton 
#'           monoculture FCM experiment
#' @param sfts character vector of gating function to test.
#' @param channels channels to be used for gating
#' @param funargs_list additional options for the chosen gating function
#' 
#' @examples
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
#' #phytoFilter specification
#' gateFunc(flowfile = reducedFlowframe(cells_nodebris),
#'               channels = c("RED.B.HLin", "YEL.B.HLin", 
#'               "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
#'               sfts = "phytoFilter", 
#'               list(ph = 0.1, proportion = 0.90)
#'               )
#' #flowClust specification
#' gateFunc(flowfile = reducedFlowframe(cells_nodebris),
#'               channels = c("RED.B.HLin", "YEL.B.HLin", 
#'               "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
#'               sfts = "flowClust", 
#'               list(K = 1:4, B = 100)
#'               )
#' #cytometree specification
#' gateFunc(flowfile = reducedFlowframe(cells_nodebris),
#'               channels = c("RED.B.HLin", "YEL.B.HLin", 
#'               "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
#'               sfts = "cytometree", 
#'               list(minleaf = 1, t = 0.10)
#'               )
#'
#' @return a flowFrame with cluster indicator generated by the software used
#'         added to the expression matrix.
#'
#' @description This function gates all flowFrames in the supplied flowSet
#'              to attach cluster labels. Then it mixes up the flowSet into
#'              one giant flowFrame and re-gates this to attach another label.
#'              These labels are used to examine if the gating algorithms can
#'              reproduce the earlier clusters before the mixing.
#' @export gateFunc

gateFunc <- function(flowfile, sfts = c("phytoFilter",
                                 "flowClust",
                                 "cytometree"), channels,
                    funargs_list
    ) {
   #columns <- c(colnames(fs[[1]]), "Clusters")
   
   nms <- names(funargs_list)
   if(sfts == 'phytoFilter') {
      
      ph <- ifelse('ph' %in% nms,
                   unlist(funargs_list[nms == 'ph']),
                   0.1
      )
      #ph <- ifelse(is.na(ph), 0.1, ph)
      proportion <-ifelse('proportion' %in% nms, 
                          unlist(funargs_list[nms == 'proportion']),
                          0.8
      )
      #proportion <- ifelse(is.na(proportion), 0.8, proportion)
      mono_gate <- phytoFilter(flowfile = flowfile, 
                               pig_channels = channels,
                               com_channels = NULL,
                               ph = ph,
                               proportion = proportion)
      f_new <- reducedFlowframe(mono_gate)
     
      
   } else if(sfts == 'flowClust') {
      
      if('K' %in% nms) {
         K <- unlist(funargs_list[nms == 'K'])
      } else {
         K <- seq_len(5)
      }
      
      
      B <- ifelse('B' %in% nms, 
                  unlist(funargs_list[nms == 'B']),
                  200
      )
      
      z_cutoff <- ifelse('z.cutoff' %in% nms, 
                  unlist(funargs_list[nms == 'z.cutoff']),
                  0.8
      )
      mono_gate <- flowClust::flowClust(flowfile, 
                                        varNames = channels,
                                        K = K, 
                                        B = B)
      min_f <- mono_gate@index
      flowClust::ruleOutliers(mono_gate[[min_f]]) <- list(z.cutoff = z_cutoff)
      
      grp_c <- apply(mono_gate[[min_f]]@z, 1, function(x) {
         if(sum(is.na(x)) < length(x)){
            return(which(x == max(x)))
         } else { 
            return(0)
         }
      })
      grp_c[is.na(grp_c)] <- 0
      f_new <- newFlowframe(flowfile, grp_c)
      
   } else if(sfts == 'cytometree') {
      
      minleaf <- ifelse('minleaf' %in% nms, 
                        unlist(funargs_list[nms == 'minleaf']),
                        1
      )
      #minleaf <- ifelse(is.na(minleaf), 1, minleaf)
      t <- ifelse('t' %in% nms, 
                  unlist(funargs_list[nms == 't']),
                  0.1
      )
      #t <- ifelse(is.na(t), 0.1, t)
      mono_gate <- cytometree::CytomeTree(exprs(flowfile)[, channels], 
                                          minleaf = , 
                                          t = t)
      f_new <- newFlowframe(flowfile, group = mono_gate$labels)
      
   } else stop('Check your input')
   return(f_new)
   
}


#' tests the accuracy of several automated gating functions on monoculture
#' flow cytometry experiments.
#'
#' @param fs flowSet with each flowFrame being a phytoplankton 
#'           monoculture FCM experiment
#' @param sfts character vector of gating function to test.
#' @param channels channels to be used for gating
#' @param nrun number of times the resampling should be done
#' @param ... extra options to be parsed to the gating function
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
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
#' #phytoFilter specification
#' gateFunc(flowfile = reducedFlowframe(cells_nodebris),
#'               channels = c("RED.B.HLin", "YEL.B.HLin", 
#'               "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
#'               sfts = "phytoFilter", 
#'               list(ph = 0.1, proportion = 0.90)
#'               )
#' @return a named list containing the following objects; \itemize{
#' \item \strong{depth -} the multivariate-depth (median) of each
#' flowFrame in the flowset supplied 
#' \item \strong{accuracy -} computed accuracy based on resampling after
#' joining the flowFrames together.
#' } 
#'
#' @description This function gates all flowFrames in the supplied flowSet
#'              to attach cluster labels. Then it mixes up the flowSet into
#'              one giant flowFrame and re-gates this to attach another label.
#'              These labels are used to examine if the gating algorithms can
#'              reproduce the earlier clusters before the mixing.
#' @export accTest


accTest <- function(fs, sfts = c("phytoFilter",
                                 "flowClust",
                                 "cytometree"), 
                    channels, 
                    nrun = 10000, ...) {
  
   funargs_list <- list(...)
   
   #gate each flowframe in the flowSet
   sft_type <- match.arg(sfts)
   
   if(is.flowSet(fs)) {
      
      dGates <- flowCore::fsApply(fs, function(x) {
         gateFunc(x, sft_type, channels, funargs_list)
      })
      
      meds <- unlist(flowCore::fsApply(fs, function(x) {
         mrfDepth::hdepthmedian(flowCore::exprs(x)[,channels])$median
      }))
      
      ## remove all expression matrix
      matss <- vector("list", length(dGates))
      for(i in seq_len(length(dGates))){
         x <- dGates[[i]]
         fexp2 <- cbind(flowCore::exprs(x), Type = rep(i, nrow(x)))
         colnames(fexp2)[which(colnames(fexp2) == "Clusters")] <- "Clusters_Mono"
         matss[[i]] <- fexp2
      }
      
      matss2 <- do.call(rbind, matss)
      
      ## form the new flowframe combining the flowFrames in the flowSet fs
      dvarMetadata <- flowCore::varMetadata(flowCore::parameters(fs[[1]]))
      ddimnames <- Biobase::dimLabels(flowCore::parameters(fs[[1]]))
      describe <- flowCore::keyword(fs[[1]])
      pd <- pData(flowCore::parameters(fs[[1]]))
      ddata <- data.frame(rbind(pd,
                          c("Clusters_Mono", "Monoculture Cluster Label", 
                            min(matss2[, 'Clusters_Mono']), 
                            max(matss2[, 'Clusters_Mono']), 
                            diff(range(max(matss2[, 'Clusters_Mono'])))
                            ),
                          c("Type", "Monoculture Type", 
                            min(matss2[, 'Type']), 
                            max(matss2[, 'Type']), 
                            diff(range(max(matss2[, 'Type'])))
                            )
                          ))
      
      paraa <- Biobase::AnnotatedDataFrame(data = ddata, 
                                           varMetadata = dvarMetadata,
                                           dimLabels = ddimnames)
      comb_fs <- flowFrame(exprs = matss2, 
                           parameters = paraa,
                           description = describe
                         )
      #gate the joined flowFrame
      mult_gates <- gateFunc(comb_fs, sft_type, channels, funargs_list)
      accs <- accuracy(mat = flowCore::exprs(mult_gates), 
                       mono_clust = "Clusters_Mono", 
                       bi_clust = "Clusters", nrun = nrun)
      accs_comp <- sum(accs == 2 | accs == 0)/ nrun
      return(list(
         depth = meds,
         accuracy = accs_comp)
         )
      
   } else stop('Object is not a flowSet')
}

#' samples two rows in a matrix and check if the samples are similar or different
#' based on their cluster labels
#'
#' @param mat matrix to be sampled from
#' @param mono_clust monoculture cluster label
#' @param bi_clust biculture cluster label
#' @param nrun number of times the resampling should be carried out.
#'             Defaults to 10000
#'
#' @examples
#' x <- matrix(NA, nrow = 100, ncol = 3)
#' xx <- apply(x, 2, rnorm, 100)
#' xx <- cbind(xx, Mono = rep(1:2, each = 50),
#'             Bi = rep(1:2, times = 50))
#' accuracy(xx, "Mono", "Bi", nrun = 5000)
#' @return a vector of integer values
#'
#' @description This function 
#' @export accuracy

accuracy <- function(mat, mono_clust, bi_clust, nrun = 10000) {
   
   tss <- numeric(nrun)
   
   dwin <- function(x) {diff(range(x)) < .Machine$double.eps ^ 0.5}
   
   for(j in seq_len(nrun)){
      
      parts <- sample(seq_len(nrow(mat)), 2, replace = FALSE)
      tests <- apply(mat[parts, c(mono_clust, bi_clust)], 2, dwin)
      tss[j] <- sum(tests)
      
   }
   #tss = 2 or tss = 0 means both the algorithm agrees for both monoculture and
   #when it is mixed.
   return(tss)
}
