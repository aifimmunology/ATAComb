#' Binarize/convert all values in a sparse matrix to 1
#'
#' @param mat a sparse matrix object from the Matrix package, e.g. dgCMatrix.
#'
#' @return a sparse matrix object
#' @export
binarize <- function(mat) {
  mat@x <- rep(1, length(mat@x))
  mat
}
