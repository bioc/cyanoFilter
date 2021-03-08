#' Removes NA values from the expression matrix of a flow cytometer file.
#'
#' @param x flowframe with expression matrix containing NAs.
#' @return flowframe with expression matrix rid of NAs.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#' mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' noNA(x = flowfile)
#'
#'
#' @importFrom methods new
#' @export noNA

noNA <- function(x) {
    dtest <- !apply(flowCore::exprs(x), 1, 
                    function(row) any(is.na(row) | is.nan(row)))
    exx <- flowCore::exprs(x)[dtest == TRUE, ]
    paraa <- flowCore::parameters(x)
    describe <- flowCore::keyword(x)
    return(flowCore::flowFrame(exprs = exx, 
                               parameters = paraa, 
                               description = describe))
}
