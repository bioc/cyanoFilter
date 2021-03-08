#' takes a list of flowframes or a flowSet with similar annotations but 
#' different expression matrices and stacks them.
#'
#' @param flowframes a list of flowframes or a flowSet
#' @return flowframe with indicators for particle cluster
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' oneDgate(flowfile, 'RED.B.HLin')
#'
#' @export newFlowframe2

newFlowframe2 <- function(flowframes) {

  mode <- match.arg(mode)

  if(methods::is(flowframes, "list")) {

    mats <- do.call(rbind, lapply(flowframes, exprs))

  } else if(methods::is(flowframes, "flowSet")) {

    ats <- do.call(rbind, fsApply(flowframes, exprs))

  } else stop("invalid object supplied")


  pd <- pData(flowCore::parameters(flowframes[[1]]))
  dvarMetadata <- flowCore::varMetadata(flowCore::parameters(flowframes[[1]]))
  ddimnames <- dimLabels(flowCore::parameters(flowframes[[1]]))
  
  paraa <- Biobase::AnnotatedDataFrame(data = pd,
                                       varMetadata = dvarMetadata,
                                       dimLabels = ddimnames)
  describe <- flowCore::keyword(flowframes[[1]])

  # full flow frame with indicator for particly type
  fflowframe <- flowCore::flowFrame(exprs = mats, parameters = paraa,
                          description = describe)
  return(fflowframe)

}
