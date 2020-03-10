#' gates out or assign indicators to phytoplankton cells based on the expression of
#' the measured pigments.
#'
#' @param flowfile flowframe after debris are removed.
#'
#' @param com_channels flowcytometer channels measuring cell complexity.
#'
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
#' }
#'
#' @description This function takes in a flowframe with debris removed and identifies
#'              the different phytoplankton cell population based on cell complexity.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the \emph{flowDensity} package to
#'          identify peaks and identify cut-off points between these peaks. This function is
#'          not designed to be called in isolation, if called in isolation an error will be
#'          returned. It is preferably called on the results from \code{\link{debris_nc}} or
#'          \code{\link{debris_inc}} function. A graph with horizontal
#'          and vertical lines used in separating the populations is returned and
#'          if \emph{to_retain = "refined"}, a circle made of dashed lines is drawn around
#'          phytoplankton cell population points.
#'
#' @seealso
#'
#'
#' @importFrom utils capture.output
#' @export pigment_gate

complexity_gate <- function(flowfile, com_channels, ph = 0.05, width_channel = "SSC.W") {

  ########### Determine which pigments needs to be gated
  comp_gates <- sapply(com_channels, get_channel, flowfile = flowfile, ph = ph)

  #channels to be gated, i.e. channels with more than one peak along that pigment
  if (sum(!is.na(comp_gates)) == 0) {

    togate_comp <- NA

  } else {

    togate_comp <- comp_gates[!is.na(comp_gates)]
  }

  if(sum(!is.na(togate_comp)) != 0 & length(togate_comp) == 1) {

    #only one channel to gate along the pigments
    cluster_gates <- oneDgate(flowfile = flowfile, togate = togate_comp,
                              width_channel = width_channel)

    #full flowframe to be returned
    full_flowframe <- new_flowframe(flowfile, group = cluster_gates$phy_ind, togate = togate_comp)

    # a vector with same length as flowframe containing the cluster each particle belongs to
    phy_ind <- cluster_gates$phy_ind

    #reduced flowframes
    red_frames <- lapply(cluster_gates$phy_red_pos, function(x) {

      full_flowframe[x, ]
    })

    reduced_frame <- as(red_frames, "flowSet")

    #construct new flowframe based on the width and togate channel
    #reduced_flowframe <- new_flowframe2(red_frames)


  } else if(sum(!is.na(togate_comp)) != 0 & length(togate_comp) > 1) {

    #more than one pigment channel to gate along the pigments
    full_flowframe <- flowfile
    # a vector with same length as flowframe containing the cluster each particle belongs to
    phy_ind <- reduced_frame <- vector("list", length = length(togate))

    for(i in 1:length(togate_comp)) {

      cluster_gates <- oneDgate(flowfile = flowfile, togate = togate_comp[i],
                                width_channel = width_channel)

      #full flowframe to be returned
      full_flowframe <- new_flowframe(flowfile = full_flowframe, group = cluster_gates$phy_ind,
                                      togate = togate_comp[i])

      #
      phy_ind[[i]] <- cluster_gates$phy_ind

      #reduced flowframes
      red_frames <- lapply(cluster_gates$phy_red_pos, function(x) {

        ff <- new_flowframe(flowfile = flowfile,
                            group = cluster_gates$phy_ind,
                            togate = togate[i])
        ff[x, ]
      })

      reduced_frame[[i]] <- as(red_frames, "flowSet")

    }
    names(phy_ind) <- names(reduced_frame) <- togate


  } else {

    #nothing to gate along the pigments
    message("Only one peak across the size and complexity channels.")
    cluster_groups <- phy_pos
    ret_flowframe <- flowfile

  }

  return(list(full_flowframe = full_flowframe,
              phy_ind = phy_ind,
              gated_channels = togate_comp,
              reduced_flowframe = reduced_frame))

}



