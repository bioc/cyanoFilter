#' cyanoFilter: A package to identify and cluster phytoplankton cells contained in flow cytometry data.
#'
#' The package provides two categories of functions:
#' \emph{metafile} preprocessing functions and \emph{fcsfile} processing functions.
#'
#' @section  metafile preprocessing functions:
#'           This set of functions (\code{\link{goodfcs}} and \code{\link{retain}}) helps to
#'           identify the appropriate fcs file to read.
#'
#' @section fcsfile processing functions:
#'          These functions (\code{\link{nona}} and \code{\link{noneg}}, \code{\link{noneg}},
#'          \code{\link{phyto_filter}}, \code{\link{celldebris_emclustering}})
#'          works on the fcs file to identify the phytoplankton populations contained in
#'          the fcs file.
#'
#' @docType package
#'
#' @name cyanoFilter

NULL
