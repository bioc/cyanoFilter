#' returns the labels stating the cluster of each row in a flowfile.
#'
#' @param flowfile flowframe after debris are removed.
#' @param togate channels detected to have more than one peak present. Provide by the
#'               \code{\link{get_channel}} function.
#'
#'
#' @export oneDgate


oneDgate <- function(flowfile, togate) {

  gates <- flowDensity::deGate(flowfile, togate, all.cuts = TRUE)

  #pop_rows is a list containing at least two vectors
  pop_rows <- row_numbers(flowframe = flowfile, gates = gates, ch = togate)

  ########## assign labels to each group
  phy_ind <- rep(NA, times = nrow(flowfile))

  for(i in 1:length(pop_rows)) {

    phy_ind[pop_rows[[i]]] <- i

  }

  return(
    list(
      phy_ind = phy_ind
        )
    )
}
