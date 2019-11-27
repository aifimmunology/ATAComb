#' Count overlaps between each GRanges object in a list and a single target GRanges object (reference).
#'
#' @param query_fragments A GenomicRanges object or a list of GenomicRanges objects to use as a query
#' @param target_GRanges A single GenomicRanges object to use as a set of target regions
#' @param binarize A logical object indicating whether or not to binarize overlaps (count any number of overlaps between query and a given target region as 1). Default is TRUE.
#' @param sparse A logical object indicating whether the results should be a sparse matrix (TRUE) or a full matrix (FALSE). Default is TRUE.
#' @param aggregate A logical object indicating whether the results should be aggregated to a vector with the sum of counts for each query to all target regions. Default is FALSE.
#' @param n_threads A numeric object specifying the number of threads to use. Default is 1.
#'
#' @return If aggregate = TRUE, a vector of counts. If aggregate is FALSE and sparse is FALSE, a matrix object. If sparse is true, as dgCMatrix object.
#' @export
count_frag_ol_ref <-function (query_fragments,
                              target_GRanges,
                              binarize = TRUE,
                              sparse = FALSE,
                              aggregate = FALSE,
                              n_threads = 1) {

  if(class(fragments) != "list") {
    fragments <- list(query_fragments = fragments)
  }

  if (aggregate) {
    out <- vector(length(fragments))
    names(out) <- names(fragments)
  } else {
    if (sparse) {
      out <- Matrix::sparseMatrix(i = integer(0),
                                  j = integer(0),
                                  dims = c(length(target_GRanges),
                                           length(fragments)))
      out <- as(out, "dgCMatrix")
    } else {
      out <- matrix(nrow = length(target_GRanges), ncol = length(fragments))
    }
    rownames(out) <- names(target_GRanges)
  }

  fragment_counts <- mclapply(1:length(fragments),
                              function(frags) {
                                ol <- GenomicRanges::countOverlaps(target_GRanges,
                                                                   fragments[[frags]])
                                frag_name <- names(fragments)[frags]

                                if(aggregate) {
                                  if(binarize) {
                                    total_count = sum(ol > 0)
                                  } else {
                                    total_count = sum(ol)
                                  }
                                  list(fragment_name = frag_name,
                                       total_count = total_count)
                                } else {
                                  if (sparse) {
                                    i <- which(ol > 0)
                                    new_count <- length(i)
                                    if(binarize) {
                                      x <- rep(1, new_count)
                                    } else {
                                      x <- ol[i]
                                    }
                                    i <- i - 1L
                                    list(fragment_name = frag_name,
                                         x = x,
                                         i = i,
                                         n_vals = new_count)
                                  } else {
                                    if(binarize) {
                                      ol[ol > 0] <- 1
                                    }
                                    list(fragment_name = frag_name,
                                         counts = ol)
                                  }
                                }
                              },
                              mc.cores = n_threads)

  if(aggregate) {
    for(fc in 1:length(fragment_counts)) {
      out[fragment_counts[[fc]]$fragment_name] <- fragment_counts[[fc]]$total_count
    }
  } else {
    frag_names <- unlist(lapply(fragment_counts,
                                function(fc) fc$fragment_name))

    if(sparse) {
      nv <- unlist(lapply(fragment_counts,
                          function(fc) fc$n_vals ))

      len_x <- sum(nv)
      out@i <- unlist(lapply(fragment_counts,
                             function(fc) fc$i))
      out@x <- as.numeric(unlist(lapply(fragment_counts,
                                        function(fc) fc$x)))
      out@p <- c(0L, as.integer(cumsum(nv)))
    } else {
      out[1:length(out)] <- unlist(lapply(fragment_counts,
                                          function(fc) fc$counts))
    }

    colnames(out) <- frag_names
    out <- out[, names(fragments)]
  }

  return(out)
}
