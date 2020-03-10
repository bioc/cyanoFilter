#' takes a flowframes, a group indicator and formulates another flowframe with group indicator as part of the
#' expression matrix of the new flowframe.
#'
#' @param flowfile flowframe after debris are removed.
#' @param group cluster group to be added to the expression matrix
#' @param togate channel detected to have more than one peak based on \code{\link{get_channel}}
#'
#' @export new_flowframe

new_flowframe <- function(flowfile, group, togate) {


  ddata <- data.frame(name = paste(togate, "Cluster", sep = "_"),
                      desc = paste("Cluster Group for", togate),
                      range = max(group) - min(group),
                      minRange = range(group)[1],
                      maxRange = range(group)[2]
                      )
  ddata <- rbind(flowfile@parameters@data, ddata)
  row.names(ddata) <- c(row.names(flowfile@parameters@data),
                        paste("$P", length(row.names(flowfile@parameters@data))+1,
                              sep = ""))


  dvarMetadata <- flowfile@parameters@varMetadata
  ddimnames <- flowfile@parameters@dimLabels
  paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                       dimLabels = ddimnames)
  describe <- flowfile@description


  #flowframe with indicator added for each cluster
  nexp_mat <- as.matrix(cbind(flowCore::exprs(flowfile), group))
  # giving a name to the newly added column to the expression matrix
  colnames(nexp_mat) <- ddata$name

  # full flow frame with indicator for particly type
  fflowframe <- methods::new("flowFrame", exprs = nexp_mat, parameters = paraa,
                             description = describe)
  return(fflowframe)
}
