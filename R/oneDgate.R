#' returns the labels stating the cluster of each row in a flowfile.
#'
#' @param flowfile flowframe after debris are removed.
#' @param togate channels detected to have more than one peak present. Provide by the
#'               \code{\link{get_channel}} function.
#' @param width_channel a second channel with only one peak present.
#'
#'
#' @export oneDgate


oneDgate <- function(flowfile, togate, width_channel) {

  gates <- flowDensity::deGate(flowfile, togate, all.cuts = TRUE)

  #pop_rows is a list containing at least two vectors
  pop_rows <- row_numbers(flowfile = flowfile, gates = gates, ch = togate)

  ########## assign labels to each group
  phy_ind <- rep(NA, times = nrow(flowfile))

  for(i in 1:length(pop_rows)) {

    phy_ind[pop_rows[[i]]] <- i

  }

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
}
