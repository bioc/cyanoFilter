#' returns the position of the identified clusters after 2D gating with the size channel
#'
#' @param flowfile flowframe after debris are removed.
#' @param phy_pos row numbers for events belonging to each cluster. This is provided by the
#'               \code{\link{row_numbers}} function.
#' @param other_channel extracted as the width_channel in \code{\link{oneDgate}}
#' @param togate channels detected to have more than one peak present. Provide by the
#'               \code{\link{get_channel}} function.
#' @param percent cut-off percentile for other_channel.
#'
#' @export refine_gate

refine_gate <- function(flowfile, phy_pos,
                        other_channel,
                        togate, percent = 0.60) {

  #solu <- vector("list", length(other_channels))


      ref_gates <- flowDensity::flowDensity(flowfile[phy_pos, ], channels = c(other_channel, togate),
                                            position = c(F, F),
                                            use.upper = c(F, T),
                                            upper = c(NA, T),
                                            use.percentile = c(T, F),
                                            percentile = c(percent, NA),
                                            ellip.gate = T)
      cynb_pos <- phy_pos[ref_gates@index]
      hulls <- ref_gates@filter

      return(list(cynb_pos = cynb_pos, hulls = hulls))

}
