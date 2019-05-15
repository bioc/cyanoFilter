#' indicates if measurement from a flowfile is good or bad @param fcsfile flowfile to be checked.
#'
#' @param metafile associated metafile to the supplied fcsfile. This is a csv file containig computed stats from the flow cytometer.
#' @param col_cpml column name or column number in metafile containing cell per microlitre measurements.
#' @param mxd_cellpML maximal accepted cell per microlitre. Flowfiles with larger cell per microlitre are termed bad.
#'                    Defaults to 1000.
#' @param mnd_cellpML minimum accepted cell per microlitre. Flowfiles with lesser cell per microlitre are termed bad.
#'                    Defaults to 50.
#'
#' @return character vector with length same as the number of rows in the metafile and its entries are
#'          \strong{good} for good files and \strong{bad} for bad files.
#'
#' @examples \dontrun{
#' goodfcs(metafile = flowframe, col_cpml = "CellspML", mxd_cellpML = 1000, mnd_cellpML = 50)
#' }
#' @export goodfcs

goodfcs <- function(metafile, col_cpml = "CellspML", mxd_cellpML = 1000, mnd_cellpML = 50) {
    if (!is.null(metafile) & !is.null(mxd_cellpML & !is.null(mnd_cellpML))) {

        goodfile <- ifelse((metafile[, col_cpml] < mxd_cellpML & metafile[, col_cpml] > mnd_cellpML), "good", "bad")

    } else stop("At least metafile is empty")
    return(goodfile)
}
