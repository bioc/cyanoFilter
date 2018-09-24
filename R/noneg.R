#' Removes negative values from the expression matrix
noneg <- function(x){
  dtest <- !apply(flowCore::exprs(x), 1 ,function (row) any(row <= 0))
  exx <- flowCore::exprs(x)[dtest == T, ]
  paraa <- x@parameters
  describe <- x@description
  methods::new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
