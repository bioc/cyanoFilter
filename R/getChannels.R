#' returns the channel with more than one peak present. 
#' It returns NA if there is only one peak present.
#'
#' @param flowfile flowframe after debris are removed.
#' @param ch channel to be checked for multiple peaks.
#' @param ph maximum peak height to be ignored. This allows 
#'           ignoring of tiny peaks that could
#'           affect the gating process.
#'
#' @return name of channel with more than one peak
#' 
#'  @examples 
#'  flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' get_channel(flowfile_logtrans, 'RED.B.HLin', 0.05) 
#'
#' @export get_channel

get_channel <- function(flowfile, ch, ph) {

  ptt <- flowDensity::getPeaks(flowfile, ch, tinypeak.removal = ph)
  if(length(ptt$Peaks) == 1) {
    return('')
  } else return(ch)

}
