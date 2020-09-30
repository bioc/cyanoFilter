
#' function to check if object is of class cyanoFilter
#' @param x any R object
#' @export is.cyanoFilter

is.cyanoFilter <- function(x) {

  if(stringr::str_detect(class(x), 'cyanoFilter')) {
    ret <- TRUE
  } else {
    ret <- FALSE
  }
  return(ret)
}



