#' Get ArchR matrices in dgCMatrix format
#'
#' @param proj an ArchRProject object
#' @param target The name of a matrix in the ArchRProject. Options are: "GeneScoreMatrix", "PeakMatrix", "TileMatrix"
#' @param ... Additional parameters passed to ArchR::getMatrixFromProject()
#'
#' @return a dgCMatrix object
#' @export
#'
#' @examples
get_archr_dgCMatrix <- function(proj,
                                target = "GeneScoreMatrix",
                                ...) {

  if(matrix == "GeneScoreMatrix") {

    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "GeneScoreMatrix",
      ...)

    mat <- se@assays@data@listData$GeneScoreMatrix

    rownames(mat) <- se@elementMetadata@listData$name

  } else if(matrix == "PeakMatrix") {
    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "PeakMatrix",
      ...)

    mat <- se@assays@data@listData$PeakMatrix

    rownames(mat) <- paste(
      seqnames(se@rowRanges),
      start(se@rowRanges),
      end(se@rowRanges),
      sep = "_")

  } else if(matrix == "TileMatrix") {
    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "TileMatrix",
      binarize = TRUE,
      ...)

    mat <- se@assays@data@listData$TileMatrix

    rownames(mat) <- paste(
      seqnames(se@rowRanges),
      start(se@rowRanges),
      end(se@rowRanges),
      sep = "_")
  }

  mat

}

