#' gates out or assign indicators to phytoplankton cells based on the expression of
#' the measured pigments.
#'
#' @param flowfile flowframe after debris are removed.
#'
#' @param pig_channels flowcytometer channels measuring phytoplankton pigmentations.
#'
<<<<<<< HEAD
#' @param ph maximum peak height to be ignored. This allows ignoring of tiny peaks that could
#'           affect the gating process.
#'
#'
#' @return list containing; \itemize{
#' \item \strong{full_flowframe -} flowframe containing only phytoplankton cells
#' \item \strong{phy_ind -} indicator for phytoplankton clusters found
#' \item \strong{gated_channels -} pigment channels with more than one peak
=======
#'
#' @param ph maximum peak height to be ignored. This allows ignoring of tiny peaks that could
#'           affect the gating process.
#'
#' @param width_channel a second channel with only one peak present.
#'
#'
#' @return list containing; \itemize{
#' \item \strong{full_flowframe -} flowframe containing only phytoplankton cells
#' \item \strong{phy_pos_nk -} unidentified particle positions
#' \item \strong{syn_pos -} phytoplankton cells positions
#' \item \strong{phy_pos_nk2 -} other unidentified particle positions
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
#' }
#'
#' @description This function takes in a flowframe with debris removed and identifies
#'              phytoplankton cell population in the provided frame.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the \emph{flowDensity} package to
<<<<<<< HEAD
#'          identify peaks and identify cut-off points between these peaks.
=======
#'          identify peaks and identify cut-off points between these peaks. This function is
#'          not designed to be called in isolation, if called in isolation an error will be
#'          returned. It is preferably called on the results from \code{\link{debris_nc}} or
#'          \code{\link{debris_inc}} function. A graph with horizontal
#'          and vertical lines used in separating the populations is returned and
#'          if \emph{to_retain = "refined"}, a circle made of dashed lines is drawn around
#'          phytoplankton cell population points.
#'
#' @seealso
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
#'
#'
#' @importFrom utils capture.output
#' @export pigment_gate

<<<<<<< HEAD
pigment_gate <- function(flowfile, pig_channels, ph = 0.05) {
=======
pigment_gate <- function(flowfile, pig_channels, ph = 0.05, width_channel = "SSC.W") {
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

  ########### Determine which pigments needs to be gated
 pot_gates <- sapply(pig_channels, get_channel, flowfile = flowfile, ph = ph)

 #channels to be gated, i.e. channels with more than one peak along that pigment
 if (sum(!is.na(pot_gates)) == 0) {

      togate <- rep(NA, length(pot_gates))

   } else {

     togate <- pot_gates[!is.na(pot_gates)]
   }

 #################### gate each channel in togate

 if(sum(!is.na(togate)) != 0 & length(togate) == 1) {

   #only one channel to gate along the pigments
<<<<<<< HEAD
   cluster_gates <- oneDgate(flowfile = flowfile, togate = togate)
=======
   cluster_gates <- oneDgate(flowfile = flowfile, togate = togate,
                              width_channel = width_channel)
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

   #full flowframe to be returned
   full_flowframe <- new_flowframe(flowfile, group = cluster_gates$phy_ind, togate = togate)

   # a vector with same length as flowframe containing the cluster each particle belongs to
   phy_ind <- cluster_gates$phy_ind

<<<<<<< HEAD
=======
   #reduced flowframes
    red_frames <- lapply(cluster_gates$phy_red_pos, function(x) {

      full_flowframe[x, ]
    })

    reduced_frame <- new_flowframe2(red_frames, mode = "list")#as(red_frames, "flowSet")

   #construct new flowframe based on the width and togate channel
   #reduced_flowframe <- new_flowframe2(red_frames)

>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

 } else if(sum(!is.na(togate)) != 0 & length(togate) > 1) {

   #more than one pigment channel to gate along the pigments
    #full flowframe to be returned
    full_flowframe <- flowfile

    # a list containing vectors with same length as flowframe containing the cluster each particle belongs to
    phy_ind <- reduced_frame <- vector("list", length = length(togate))

    for(i in 1:length(togate)) {

<<<<<<< HEAD
       cluster_gates <- oneDgate(flowfile = flowfile, togate = togate[i])
=======
       cluster_gates <- oneDgate(flowfile = flowfile, togate = togate[i],
                                 width_channel = width_channel)
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
       #full flowframe to be returned
       full_flowframe <- new_flowframe(flowfile = full_flowframe,
                                       group = cluster_gates$phy_ind,
                                       togate = togate[i])
       #cluster group indicator
       phy_ind[[i]] <- cluster_gates$phy_ind

<<<<<<< HEAD

    }
    names(phy_ind) <- togate
=======
       #reduced flowframes
       red_frames <- lapply(cluster_gates$phy_red_pos, function(x) {

          ff <- new_flowframe(flowfile = flowfile,
                              group = cluster_gates$phy_ind,
                              togate = togate[i])
          ff[x, ]
       })

       reduced_frame[[i]] <- new_flowframe2(red_frames, mode = "list")#as(red_frames, "flowSet")


    }
    names(phy_ind) <- names(reduced_frame) <- togate
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a


 } else {

   #nothing to gate along the pigments
<<<<<<< HEAD
   message("Only one peak across all the supplied channels.")
   #cluster_groups <- phy_pos
    phy_ind <- NA
    full_flowframe <- flowfile
    #reduced_frame <- NA
=======
   message("Only one peak across all supplied pigment channels.")
   #cluster_groups <- phy_pos
    full_flowframe <- flowfile
   phy_ind <- rep(NA, nrow(flowfile))
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

 }

 return(list(full_flowframe = full_flowframe,
             phy_ind = phy_ind,
<<<<<<< HEAD
             gated_channels = togate
             )
        )
=======
             gated_channels = togate,
             reduced_flowframe = reduced_frame))
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

}



