#' gates out or assign indicators to phytoplankton cells based 
#' on the expression of
#' the measured pigments.
#'
#' @param flowfile flowframe after debris are removed.
#'
#' @param pig_channels flowcytometer channels measuring phytoplankton 
#'                     pigmentations.
#'
#' @param ph maximum peak height to be ignored. This allows ignoring of 
#'           tiny peaks that could
#'           affect the gating process.
#'
#'
#' @return list containing; \itemize{
#' \item \strong{full_flowframe -} flowframe containing only phytoplankton cells
#' \item \strong{phy_ind -} indicator for phytoplankton clusters found
#' \item \strong{gated_channels -} pigment channels with more than one peak
#' }
#'
#' @description This function takes in a flowframe with debris removed and 
#'              identifies
#'              phytoplankton cell population in the provided frame.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the 
#'          \emph{flowDensity}  package to identify peaks and identify cut-off 
#'          points between these peaks.
#'          
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
#' pigment_gate(flowfile = flowfile_logtrans,
#' pig_channels = c("RED.B.HLin", "YEL.B.HLin",
#'                     "FSC.HLin", "RED.R.HLin"),
#' ph = 0.06)
#'
#'
#' @importFrom utils capture.output
#' @export pigment_gate

pigment_gate <- function(flowfile, pig_channels, ph = 0.05) {

  ########### Determine which pigments needs to be gated
 pot_gates <- vapply(pig_channels, get_channel, 
                     c(''), 
                     flowfile = flowfile, ph = ph)
 #turn empty strings to NA
 pot_gates[pot_gates == ''] <- NA

 #channels to be gated, i.e. channels with more than one peak along that pigment
 if (sum(!is.na(pot_gates)) == 0) {

      togate <- rep(NA, length(pot_gates))

   } else {

     togate <- pot_gates[!is.na(pot_gates)]
   }

 #################### gate each channel in togate

 if(sum(!is.na(togate)) != 0 & length(togate) == 1) {

   #only one channel to gate along the pigment
   cluster_gates <- oneDgate(flowfile = flowfile, togate = togate)
   #full flowframe to be returned
   full_flowframe <- new_flowframe(flowfile, group = cluster_gates$phy_ind, 
                                   togate = togate)

   # a vector with same length as flowframe containing the cluster each 
   # particle belongs to
   phy_ind <- cluster_gates$phy_ind


 } else if(sum(!is.na(togate)) != 0 & length(togate) > 1) {

   #more than one pigment channel to gate along the pigments
    #full flowframe to be returned
    full_flowframe <- flowfile

    # a list containing vectors with same length as flowframe containing the 
    # cluster each particle belongs to
    phy_ind <- reduced_frame <- vector("list", length = length(togate))

    for(i in seq_len(length(togate))) {

       cluster_gates <- oneDgate(flowfile = flowfile, togate = togate[i])
       #full flowframe to be returned
       full_flowframe <- new_flowframe(flowfile = full_flowframe,
                                       group = cluster_gates$phy_ind,
                                       togate = togate[i])
       #cluster group indicator
       phy_ind[[i]] <- cluster_gates$phy_ind

    }
    names(phy_ind) <- togate

 } else {

   #nothing to gate along the pigments
   message("Only one peak across all the supplied channels.")
   #cluster_groups <- phy_pos
    phy_ind <- NA
    full_flowframe <- flowfile
    #reduced_frame <- NA

 }

 return(list(full_flowframe = full_flowframe,
             phy_ind = phy_ind,
             gated_channels = togate
             )
        )

}



