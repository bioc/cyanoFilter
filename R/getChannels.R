#' returns the channel with more than one peak present. It returns NA if there is only one peak present.
#'
#' @param flowfile flowframe after debris are removed.
#' @param ch channel to be checked for multiple peaks.
#' @param ph maximum peak height to be ignored. This allows ignoring of tiny peaks that could
#'           affect the gating process.
#'
#' @export get_channel

get_channel <- function(flowfile, ch, ph) {

  ptt <- flowDensity::getPeaks(flowfile, ch, tinypeak.removal = ph)
  if(length(ptt$Peaks) != 1) {
    return(ch)
  } else return(NA)

}
