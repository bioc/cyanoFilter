# indicates if measurement from a flowfile is good or bad
#
# @param fcsfile flowfile to be checked.
# @param metafile associated metafile to the supplied fcsfile.
# @param CpermL column name in metafile containing cell per microlitre.
# @param rowinmetafile row containing measurements associated to fcsfile.
# @param d_cellnum desired minimum number of particle counts. Flowfiles with smaller particle counts are
#                  termed bad. Defaults to 90.
# @param d_cellpML maximal accepted cell per microlitre. Flowfiles with larger cell per microlitre are
#        termed bad. Defaults to 1000
# @return dataframe with columns
# \itemize{
#         \item \strong{CML -} cell per microlitre as read from the associated metafile.
#         \item \strong{status -} flowfile status, either "good" or "bad".
#         \item \strong{ID -} is the row number in meta file
#         \item \strong{Size -} is the number of cells in the FCS file.
#         }
# @examples
# \dontrun{
#   goodfcs(fcsfile = flowframe, metadtat, Cell_per_microlitre, 12, d_cellnum = 500, 1500)
# }
#
# @export goodfcs

goodfcs <- function(metafile, mxd_cellpML = 1000, mnd_cellpML = 50) {
  if(!is.null(metafile) & !is.null(mxd_cellpML & !is.null(mnd_cellpML))) {

    goodfile <- ifelse((metafile$CellspML < mxd_cellpML & metafile$CellspML > mnd_cellpML),
                       "good", "bad")

  } else stop("At least metafile is empty")
  return(goodfile)
}
