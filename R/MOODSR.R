#' MOODSR: R interface for MOODS
#'
#' @import Matrix
#' @import Biostrings
#' @import methods
#' @importFrom Rcpp sourceCpp
#' @useDynLib MOODSR
#' @docType package
#' @name MOODSR
NULL
#> NULL


.onUnload <- function (libpath) {
  library.dynam.unload("mypackage", libpath)
}
