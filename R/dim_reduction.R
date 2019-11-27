#' Perform Latent Semantic Indexing of a binary matrix
#'
#' @param atac_matrix A binary, dgCMatrix of overlaps between cells (columns) and peaks (rows).
#' @param site_frequency_threshold A minimum threshold for the fraction of cells for which each peak must be positive to be retained.
#' @param seed A numeric value to use as a random seed. Default is 3030.
#'
#' @return A VD matrix from singular value decomposition.
#' @export
atac_lsi <- function(atac_matrix,
                     site_frequency_threshold = 0.03,
                     seed = 3030) {

  num_cells_ncounted <- Matrix::rowSums(atac_matrix)
  threshold <- ncol(atac_matrix) * site_frequency_threshold

  ncounts <- atac_matrix[num_cells_ncounted >= threshold,]

  ## Normalize the data with TF-IDF
  nfreqs <- t(t(ncounts) / Matrix::colSums(ncounts))
  nfreqs@x <- log1p(nfreqs@x * 1e5)

  tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))

  tf_idf_counts@x[is.na(tf_idf_counts@x)] <- 0

  ## Do SVD
  set.seed(0)
  SVD <- irlba::irlba(tf_idf_counts, 50, 50, maxit=1000)
  d_diag <- matrix(0, nrow=length(SVD$d), ncol=length(SVD$d))
  diag(d_diag) <- SVD$d
  SVD_vd <- t(d_diag %*% t(SVD$v))
  rownames(SVD_vd) <- colnames(atac_matrix)
  colnames(SVD_vd) <- paste0('pca_', 1:ncol(SVD_vd))

  return(SVD_vd)
}

#' Compute Jaccard distances
#'
#' @param m a sparse matrix of KNN
#'
#' @return a sparse matrix of Jaccard distances.
#' @export
jaccard <-  function(m) {

  # common values:
  A <- Matrix::tcrossprod(m)
  A <- as(A, "dgTMatrix")

  # counts for each row
  b <- Matrix::rowSums(m)

  # Jacard formula: #common / (#i + #j - #common)
  x <- A@x / (b[A@i+1] + b[A@j+1] - A@x)

  A@x <- x

  return(A)
}
