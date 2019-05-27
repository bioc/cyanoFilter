#' plots the expression matrix of a flowframe. Note that, it takes some time to display the plot.
#'
#' @param flowfile flowframe to be plotted
#' @param notToPlot column in expression matrix not to be plotted
#'
#' @examples \dontrun{ pair_plot(x = flowfile, notToPlot = "TIME") }
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline points panel.smooth pairs smoothScatter text
#' @export pair_plot

pair_plot <- function(flowfile, notToPlot = c("TIME")) {
    toplot <- setdiff(flowCore::colnames(flowfile), notToPlot)
    col.palette <- colorRampPalette(c("white", "blue", "cyan", "green", "orange", "red"), space = "Lab")
    pairs(flowCore::exprs(flowfile)[, toplot], pch = ".",
          panel = function(...) smoothScatter(..., nrpoints = 0,
                  colramp = col.palette,
                  add = TRUE), gap = 0.2, main = flowCore::identifier(flowfile))
}
