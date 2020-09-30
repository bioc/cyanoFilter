#' returns the labels stating the cluster of each row in a flowfile.
#'
#' @param flowfile flowframe after debris are removed.
#' @param togate channels detected to have more than one peak present. Provide by the
#'               \code{\link{get_channel}} function.
<<<<<<< HEAD
=======
#' @param width_channel a second channel with only one peak present.
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
#'
#'
#' @export oneDgate


<<<<<<< HEAD
oneDgate <- function(flowfile, togate) {
=======
oneDgate <- function(flowfile, togate, width_channel) {
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

  gates <- flowDensity::deGate(flowfile, togate, all.cuts = TRUE)

  #pop_rows is a list containing at least two vectors
<<<<<<< HEAD
  pop_rows <- row_numbers(flowframe = flowfile, gates = gates, ch = togate)
=======
  pop_rows <- row_numbers(flowfile = flowfile, gates = gates, ch = togate)
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

  ########## assign labels to each group
  phy_ind <- rep(NA, times = nrow(flowfile))

  for(i in 1:length(pop_rows)) {

    phy_ind[pop_rows[[i]]] <- i

  }

<<<<<<< HEAD
  return(
    list(
      phy_ind = phy_ind
        )
    )
=======
  # refine the potentials cells
  cyb_refine <- lapply(pop_rows, function(x) {

    trf <- refine_gate(flowfile, phy_pos = x,
                other_channel = width_channel,
                togate = togate)
    return(trf)

  })

  # create data to plot the convex hulls
  hulls <- lapply(cyb_refine, "[[", "hulls")
  ns <- sapply(hulls, nrow)
  hulls_data <- do.call(rbind.data.frame, hulls)
  hull_data <- data.frame(hulls_data, Cluster.Group = rep(1:max(phy_ind), times = ns))

  # create list of phytoplankton position per population group
  phy_red_pos <- lapply(cyb_refine, "[[", "cynb_pos")


  return(list(phy_ind = phy_ind, hull_data = hull_data, phy_red_pos = phy_red_pos))
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
}
