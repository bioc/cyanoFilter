lnTrans <- function(x, notToTransform = c("SSC.W","TIME")){
  exx <- cbind(log(flowCore::exprs(x)[, which(!(colnames(x) %in% notToTransform ))]),
               flowCore::exprs(x)[, which(colnames(x) %in% notToTransform)])
  colnames(exx) <- colnames(x)
  paraa <- x@parameters
  describe <- x@description
  methods::new("flowFrame", exprs = exx, parameters = paraa, description = describe)
}
