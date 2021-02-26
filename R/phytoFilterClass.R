#' the phytofilter class
#' 
#' @slot fullflowframe object of class "flowFrame" same as the input flowFrame
#' @slot flowframe_proportion object of class "flowFrame" a partial flowframe 
#'                            containing a proportion of the measured particles
#' @slot clusters_proportion object of class "numeric" representing the 
#'                           proportion of particles in each cluster
#' @slot particles_per_cluster object of class "data.frame" representing the 
#'                             number of particles in each cluster
#' @slot Cluster_ind object of class "integer" representing the labels for each
#'                   cluster
#' @slot gated_channels object of class "character" representing the names of
#'                      channels with multiple peaks
#' @slot channels object of class "character" representing the names of the 
#'                channels 
#' @import methods
#' @export PhytoFilter            

PhytoFilter <- setClass('PhytoFilter',
         slots = list(
         fullflowframe = 'flowFrame',
         flowframe_proportion = 'flowFrame',
         clusters_proportion = 'numeric',
         particles_per_cluster = 'data.frame',
         Cluster_ind = 'integer',
         gated_channels = 'character',
         channels = 'character'
        )
)

#' constructor for the PhytoFilter class
#' 
#' @param fullflowframe same as the input flowFrame
#' @param flowframe_proportion a partial flowframe containing a proportion of 
#'                             the measured particles
#' @param clusters_proportion is the proportion of particles in each cluster
#' @param particles_per_cluster number of particles in each cluster
#' @param Cluster_ind labels for each cluster
#' @param gated_channels the names of channels with multiple peaks 
#' @param channels the names of all channels supplied to the function
#' @return object of class PhytoFilter
#' 
#'@examples
#'  flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = reducedFlowframe(cells_nonmargin),
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' phyto_filter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#' @export PhytoFilter
PhytoFilter <- function(fullflowframe, flowframe_proportion, 
                        clusters_proportion, particles_per_cluster,
                        Cluster_ind, gated_channels, channels) {
  
  new('PhytoFilter', fullflowframe = fullflowframe, 
      flowframe_proportion = flowframe_proportion, 
      clusters_proportion = clusters_proportion, 
      particles_per_cluster = particles_per_cluster,
      Cluster_ind = Cluster_ind, 
      gated_channels = gated_channels,
      channels = channels
    )
  
}

#' generic function for extracting the full flowframe
#' 
#' @param x an object of either class PhytoFilter, MarginEvents or DebrisFilter
#' @return generic to extract fullFlowframe

setGeneric("fullFlowframe", function(x){
  standardGeneric("fullFlowframe")
})

#' accesor method for reduced flowframe(PhytoFilter class)
#' @param x an object of class PhytoFilter
#' @return fullFlowframe method for PhytoFilter
#' @examples 
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = reducedFlowframe(cells_nonmargin),
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' phy1 <- phyto_filter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#' fullFlowframe(phy1)
#' @export

setMethod("fullFlowframe", "PhytoFilter", 
          function(x) { x@fullflowframe })

#' generic function for extracting the reduced flowframe
#' 
#' @param x an object of either class PhytoFilter, MarginEvents or DebrisFilter
#' @return generic for reduced flowFrame

setGeneric("reducedFlowframe", function(x){
  standardGeneric("reducedFlowframe")
})

#' accesor method for reduced flowframe(PhytoFilter class)
#' @param x an object of class PhytoFilter
#' @return reduced flowFrame method for PhytoFilter
#' #' @examples 
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = reducedFlowframe(cells_nonmargin),
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' phy1 <- phyto_filter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#' reducedFlowframe(phy1)
#' @export

setMethod("reducedFlowframe", "PhytoFilter", 
          function(x) { x@flowframe_proportion })

#' plot method for PhytoFilter objects
#' 
#' @param x an object of class PhytoFilter
#' @return object of calss ggplot
#' @export
setMethod("plot", "PhytoFilter", 
          function(x) {
            ggpairsDens(reducedFlowframe(x), channels = x@channels, 
                        group = "Clusters") 
})

#' takes a flowframes, a vector of channels, cluster indicator and return 
#' desired summaries per cluster
#'
#' @param object An object of class cyanoFilter to be summarised.
#' @param channels channels whose summaries are to be computed
#' @param cluster_var column name in expression matrix containing the cluter 
#'                    indicators
#' @param summary summary statistic of interest. Only mean and 
#'                variance-covariance matrix supported at the moment.
#' @param ... other arguments. Not used at the moment
#' @return list containing computed summaires
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'               package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#' Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = reducedFlowframe(cells_nonmargin),
#'                            ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' fin <- phyto_filter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' summary(object = reducedFlowframe(fin),
#'         channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'         cluster_var = "Clusters",
#'         summary = 'mean')
#'
#' @importFrom stats aggregate median sd
#' @export
setMethod("summary", "PhytoFilter", 
          function(object, channels = NULL,
                   cluster_var = "Clusters",
                   summary = c("mean", "median" , "cov", "n"), ...) {
            if(sum(summary %in% c("mean", "median", "cov", "n")) == 0) {
              
              stop("wrong summary option supplied")
              
            } else if(is.PhytoFilter(object) & (is.null(channels) | 
                                                is.null(cluster_var)) ) {
              
              stop("You must supply channels and cluster_var for objects 
                    of class phytoFilter")
            }
            
            if(is.PhytoFilter(object)) {
              
              ff <- reducedFlowframe(object)
              
              sm <- match.arg(summary)
              
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
                                 length)[, seq_len(2)]
                names(mns) <- c(cluster_var, "Number_of_particles")
                return(mns)
                
              } else if(sm == "cov") {
                
                #calculate variance per cluster
                uqs <- unique(ff@exprs[, cluster_var])
                vars_list <- vapply(uqs, function(i) {
                  
                  cov(ff@exprs[ff@exprs[, cluster_var] == i, channels])
                  
                }, outer(seq_len(2), seq_len(2)))
                attr(vars_list, 'dimnames') <- list(NULL, NULL, 
                                                    paste("Clusters", uqs, 
                                                          sep = "_"))
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
                vars_list <- vapply(uqs, function(i) {
                  
                  cov(ff@exprs[ff@exprs[, cluster_var] == i, channels])
                  
                }, outer(seq_len(2), seq_len(2)))
                attr(vars_list, 'dimnames') <- list(NULL, NULL, 
                                                    paste("Clusters", uqs, 
                                                          sep = "_"))
                names(vars_list) <- paste("Clusters", uqs, sep = "_")
                return(list(means = mns, medians = mns2, 
                            variances = vars_list))
              }
              
            } else stop('object of wrong class supplied')
            
            
          })



#' the marginEvent class
#' 
#' @slot fullflowframe object of class "flowFrame" same as the input flowFrame
#' @slot reducedflowframe object of class "flowFrame" a partial flowframe 
#'                            containing a proportion of the measured particles
#' @slot N_margin object of class "numeric" representing the 
#'                           proportion of particles in each cluster
#' @slot N_nonmargin object of class "integer" representing the 
#'                             number of particles in each cluster
#' @slot N_particle object of class "integer" representing the labels for each
#'                   cluster
#' @slot Channel object of class character representing channel measuring cell 
#'               width
#' @slot y_toplot object of class character representing plot variable 
#' @slot cut object of class numberic representing estimated inflection point or
#'           supplied cut-off point 
#'           
#' @export MarginEvents

MarginEvents <- setClass('MarginEvents',
         slots = list(
           fullflowframe = 'flowFrame',
           reducedflowframe = 'flowFrame',
           N_margin = 'numeric',
           N_nonmargin = 'integer',
           N_particle = 'integer',
           Channel = 'character',
           y_toplot = 'character',
           cut = 'numeric'
         )
)

#' constructor for the MarginEvents class
#' 
#' @param fullflowframe same as the input flowFrame
#' @param reducedflowframe a partial flowframe containing non-margin events
#' @param N_margin number of margin particles measured
#' @param N_nonmargin number of non-margine particles
#' @param N_particle total number of particles measured
#' @param Channel channel measuring the width of the particles
#' @param y_toplot another channel to use in a bivariate plot
#' @param cut the cut-off point estimated or supplied.
#' @return object of class MarginEvents
#' 
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                  package = "cyanoFilter",
#'                  mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1)
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cellmargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' 
#' 
#' @export MarginEvents

MarginEvents <- function(fullflowframe, reducedflowframe, 
                         N_margin, N_nonmargin,
                         N_particle, Channel, y_toplot, cut) {
  
  new('MarginEvents', fullflowframe = fullflowframe, 
      reducedflowframe = reducedflowframe, 
      N_margin = N_margin, 
      N_nonmargin = N_nonmargin,
      N_particle = N_particle,
      Channel = Channel,
      y_toplot = y_toplot,
      cut = cut
  )
  
}

#' accesor method for the fullflowframe (MarginEvent class)
#' @param x an object of class MarginEvents
#' @return full Flowframe method for MarginEvents
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#' Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' fullFlowframe(cells_nonmargin)
#' @export

setMethod("fullFlowframe", "MarginEvents", 
          function(x) { x@fullflowframe })

#' accesor method for reduced flowframe (MarginEvent class)
#' @param x an object of class MarginEvents
#' @return reduced Flowframe method for MarginEvents
#' @examples 
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#' Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' reducedFlowframe(cells_nonmargin)
#' @export

setMethod("reducedFlowframe", "MarginEvents", 
          function(x) { x@reducedflowframe })

#' plot method for MarginEvents objects
#' 
#' @param x an object of class MarginEvents
#' @return object of class ggplot
#' @export
setMethod("plot", "MarginEvents", 
          function(x) {
            ggplotDens(reducedFlowframe(x), c(x@Channel, x@y_toplot)) +
              geom_vline(xintercept = x@cut, linetype = 'dashed')
          })

#' takes a flowframes, a vector of channels, cluster indicator and return 
#' desired summaries per cluster
#'
#' @param object An object of class MarginEvents to be summarised.
#' @param channels channels whose summaries are to be computed
#' @param ... other arguments. Not used at the moment
#' @return list containing the required summaries
#' @examples 
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'               package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, 
#' Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' summary(object = reducedFlowframe(cells_nonmargin),
#'         channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"))
#' 
#' @export
setMethod("summary", "MarginEvents", 
          function(object, channels = NULL, ...) {
            if(is.MarginEvents(object) & is.null(channels) ) {
              
              stop("You must supply channels for objects 
                    of class MarginEvents")
            }
            
            if(is.MarginEvents(object)) {
              ff <- reducedFlowframe(object) 
              #calculate overall mean
              mns <- apply(ff@exprs[, channels],2,
                           mean,
                           na.rm = TRUE)
              #calculate overall meadian
              mns2 <- apply(ff@exprs[, channels],2,
                            median,  
                            na.rm = TRUE)
              #calculate overall variance covariance
              vars_list <- cov(ff@exprs[, channels])
              
              #Number of margin events
              nMargin <- object@N_margin
              #Number of nonNamrgin events
              nonMargin <- object@N_nonmargin
              
              return(list(means = mns, medians = mns2, 
                          variances = vars_list,
                          marginCount = nMargin,
                          nonMarginCount = nonMargin))
              
              
            } else stop('object of wrong class supplied')
          })



#' the Debris class
#' 
#' @slot fullflowframe object of class "flowFrame" same as the input flowFrame
#' @slot reducedflowframe object of class "flowFrame" a partial flowframe 
#'                            containing a proportion of the measured particles
#' @slot deb_pos object of class "numeric" representing the 
#'                           proportion of particles in each cluster
#' @slot syn_all_pos object of class "numeric" representing the 
#'                             number of particles in each cluster
#' @slot deb_cut object of class "numeric" representing the inflection point
#'               between debris and good cells.
#' @slot ch_chlorophyll objet of class "character" representing the chlorophyll
#'                      channel.
#' @slot ch_p2 object of class character to 
#' @export DebrisFilter

DebrisFilter <- setClass('DebrisFilter',
         slots = list(
           fullflowframe = "flowFrame",
           reducedflowframe = 'flowFrame',
           deb_pos = "numeric",
           syn_all_pos = 'numeric',
           deb_cut = 'numeric',
           ch_chlorophyll = 'character', 
           ch_p2 = 'character'
         )
)

#' constructor for the DebrisFilter class
#' 
#' @param fullflowframe same as the input flowFrame
#' @param reducedflowframe a partial flowframe containing non-margin events
#' @param deb_pos number of margin particles measured
#' @param syn_all_pos number of non-margine particles
#' @param deb_cut estimated inflection point between debris and good cells
#' @param ch_chlorophyll channel estimating chlorophyll level
#' @param ch_p2 plotting channel
#' @return object of class DebrisFilter
#' @export DebrisFilter

DebrisFilter <- function(fullflowframe, 
                         reducedflowframe, 
                         deb_pos, syn_all_pos, 
                         deb_cut, ch_chlorophyll, 
                         ch_p2) {
  
  new('DebrisFilter', 
      fullflowframe = fullflowframe, 
      reducedflowframe = reducedflowframe, 
      deb_pos = deb_pos, 
      syn_all_pos = syn_all_pos,
      deb_cut = deb_cut,
      ch_chlorophyll = ch_chlorophyll, 
      ch_p2 = ch_p2
  )
  
}

#' accesor method for reduced flowframe (DebrisFilter class)
#' @param x an object of class DebrisFilter
#' @return full flowFrame method for DebrisFilter
#' @export

setMethod("fullFlowframe",  "DebrisFilter", 
          function(x) { x@fullflowframe })

#' accesor method for reduced flowframe (DebrisFilter class)
#' @param x an object of class DebrisFilter
#' @return reduced flowFrame method for DebrisFilter
#' @export

setMethod("reducedFlowframe",  "DebrisFilter", 
          function(x) { x@reducedflowframe })

#' plot method for DebrisFilter objects
#' 
#' @param x an object of class DebrisFilter
#' @return object of class ggplot
#' @export
setMethod("plot", "DebrisFilter", 
          function(x) {
            y <- fullFlowframe(x)
            ggplotDens(y, c(x@ch_chlorophyll, x@ch_p2)) +
              geom_vline(xintercept = x@deb_cut, color = "red", 
                         linetype = "dashed") +
              geom_text(aes(x = mean(y@exprs[which(
                flowCore::exprs(y)[, x@ch_chlorophyll] <= x@deb_cut),
                x@ch_chlorophyll]),
                y = mean(y@exprs[
                  which(flowCore::exprs(y)[, x@ch_p2] <=
                          x@deb_cut), x@ch_p2])), inherit.aes = FALSE,
                label = paste0("Debris"), colour = "blue", size = 6)
          })

#' takes a flowframes, a vector of channels, cluster indicator and return 
#' desired summaries per cluster
#'
#' @param object An object of class DebrisFilter to be summarised.
#' @param channels channels whose summaries are to be computed
#' @param ... other arguments. Not used at the moment
#' @return list containing required summaries
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'               package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nodeb <- debris_nc(flowframe = flowfile_logtrans, 
#'           ch_chlorophyll = "RED.B.HLin",
#'           ch_p2 = "YEL.B.HLin",
#'           ph = 0.05)
#' 
#' summary(object = reducedFlowframe(cells_nodeb),
#'         channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"))
#' 
#' @export
setMethod("summary", "DebrisFilter", 
          function(object, channels = NULL, ...) {
            if(is.DebrisFilter(object) & is.null(channels)  ) {
              
              stop("You must supply channels and cluster_var for objects 
                    of class DebrisFilter")
            }
            
            if(is.DebrisFilter(object)) {
              
              ff <- reducedFlowframe(object) 
              
              #calculate overall mean
              mns <- apply(ff@exprs[, channels],2,
                           mean, na.rm = TRUE)
              #calculate overall mean
              mns2 <- apply(ff@exprs[, channels],2,
                            median,
                            na.rm = TRUE)
              #calculate overall variance covariance
              vars_list <- cov(ff@exprs[, channels])
              
              return(list(means = mns, medians = mns2, 
                          variances = vars_list))
              
              
            } else stop('object of wrong class supplied')
          })


#' function to check if object is of class cyanoFilter(PhytoFilter)
#' @param x any R object
#' @return TRUE if object is of class PhytoFilter. FALSE otherwise
#' 
#' @examples
#'  x <- c(1, 5, 4)
#'  is.PhytoFilter(x)
#' 
#' @export is.PhytoFilter

is.PhytoFilter <- function(x) {
  methods::is(x, 'PhytoFilter')
}


#' function to check if object is of class cyanoFilter(MarginEvents)
#' @param x any R object
#' @return TRUE if object is of class MarginEvents. FALSE otherwise
#' 
#' @examples
#'  x <- c(1, 5, 4)
#'  is.MarginEvents(x)
#' @export is.MarginEvents

is.MarginEvents <- function(x) {
  methods::is(x, 'MarginEvents')
}


#' function to check if object is of class cyanoFilter(DebrisFilter)
#' @param x any R object
#' @return TRUE if object is of class DebrisFilter. FALSE otherwise
#' 
#' @examples
#'  x <- c(1, 5, 4)
#'  is.DebrisFilter(x)
#' @export is.DebrisFilter

is.DebrisFilter <- function(x) {
  methods::is(x, 'DebrisFilter')
}









