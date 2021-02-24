#' returns the labels stating the cluster of each row in a flowfile.
#'
#' @param flowfile flowframe after debris are removed.
#' @param togate channels detected to have more than one peak present. 
#' Provide by the \code{\link{get_channel}} function.
#' @return list of indicators for cells above and below an estimated threshold
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
#' oneDgate(flowfile, 'RED.B.HLin')
#'
#' @export oneDgate

oneDgate <- function(flowfile, togate) {

  gates <- flowDensity::deGate(flowfile, togate, all.cuts = TRUE)

  #pop_rows is a list containing at least two vectors
  pop_rows <- row_numbers(flowframe = flowfile, gates = gates, ch = togate)

  ########## assign labels to each group
  phy_ind <- rep(NA, times = nrow(flowfile))

  for(i in seq_len(length(pop_rows))) {

    phy_ind[pop_rows[[i]]] <- i

  }

  return(
    list(
      phy_ind = phy_ind
        )
    )

}
