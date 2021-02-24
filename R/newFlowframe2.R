#' takes a list of flowframes or a flowSet with similar annotations but 
#' different expression matrices and stacks them.
#'
#' @param flowframes a list of flowframes or a flowSet
#' @param mode the type of object flowframe is. Can either be "list" or "set".
#' @return flowframe with indicators for particle cluster
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
#' oneDgate(flowfile, 'RED.B.HLin')
#'
#' @export new_flowframe2

new_flowframe2 <- function(flowframes, mode = c("list", "set")) {

  mode <- match.arg(mode)

  if(mode == "list") {

    mats <- do.call(rbind,lapply(flowframes, exprs))

  } else if(mode == "set") {

    ats <- do.call(rbind, fsApply(flowframes, exprs))

  } else stop("invalid mode supplied")


  dvarMetadata <- flowframes[[1]]@parameters@varMetadata
  ddimnames <- flowframes[[1]]@parameters@dimLabels
  paraa <- paraa <- Biobase::AnnotatedDataFrame(data = 
                                              flowframes[[1]]@parameters@data,
                                                varMetadata = dvarMetadata,
                                                dimLabels = ddimnames)
  describe <- flowframes[[1]]@description

  # full flow frame with indicator for particly type
  fflowframe <- flowCore::flowFrame(exprs = mats, parameters = paraa,
                          description = describe)
  return(fflowframe)

}
