#' Removes NA values from the expression matrix of a flow cytometer file.
#'
#' @param x flowframe with expression matrix containing NAs.
#' @return flowframe with expression matrix rid of NAs.
#'
#' @examples
#' \dontrun{
#' nona(x = flowfile)
#' }
#'
#' @importFrom methods new
#' @export nona

nona <- function(x) {
    dtest <- !apply(flowCore::exprs(x), 1, function(row) any(is.na(row) | is.nan(row)))
    exx <- flowCore::exprs(x)[dtest == T, ]
    paraa <- x@parameters
    describe <- x@description
    methods::new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
