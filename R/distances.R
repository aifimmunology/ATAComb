#' Compute Jaccard distances
#'
#' @param m a sparse matrix of KNN
#'
#' @return a sparse matrix of Jaccard distances.
#' @export
sparse_jaccard <-  function(m) {

  # common values:
  A <- Matrix::tcrossprod(m)
  A <- as(A, "dgTMatrix")

  # counts for each row
  b <- Matrix::rowSums(m)

  # Jacard formula: #common / (#i + #j - #common)
  # i.e. intersect / union
  A@x <- A@x / (b[A@i+1] + b[A@j+1] - A@x)

  return(A)
}

#' Compute Jaccard distances between a vector (i.e. a single row) against all rows of a sparse matrix
#'
#' Note: for efficiency, if axis = "row", will transpose and split the target_mat,
#' which will require available RAM equal to twice the size of target_mat.
#'
#' If axis = "col", will only split, which requires available RAM equal to the size of target_mat.
#'
#' @param query_vec A binary vector object to use as a query.
#' @param target_mat A dgCMatrix to use as targets (will be treated as binary).
#' @param axis Either "row" or "col".
#'
#' @return a vector of Jaccard distances between the query_vec and each row of the target_mat
#' @export
feature_jaccard_distance <- function(query_vec, target_mat, axis = "row") {
  #assertthat::assert_that(sum(class(query_vec) %in% c("numeric","integer") > 0))
  #assertthat::assert_that("dgCMatrix" %in% class(target_mat))
  #assertthat::assert_that(class(axis) == "character")
  #assertthat::assert_that(length(axis) == 1)

  match.arg(axis, c("col","row"))

  a <- which(query_vec == 1)

  if(axis == "row") {
    # transpose so genes/peaks are columns
    target_mat <- t(target_mat)
  }

  il <- split(target_mat@i, rep(1:ncol(target_mat), diff(target_mat@p)))
  rm(target_mat)

  unlist(
    lapply(
      il,
      function(b) {
        1 - length(intersect(b+1,a)) / length(union(b+1,a))
      }
    )
  )
}
