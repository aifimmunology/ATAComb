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

  if(class(query_fragments) != "list") {
    query_fragments <- list(query_fragments = query_fragments)
  }

  if (aggregate) {
    out <- integer(length(query_fragments))
    names(out) <- names(query_fragments)
  } else {
    if (sparse) {
      out <- Matrix::sparseMatrix(i = integer(0),
                                  j = integer(0),
                                  dims = c(length(target_GRanges),
                                           length(query_fragments)))
      out <- as(out, "dgCMatrix")
    } else {
      out <- matrix(nrow = length(target_GRanges), ncol = length(query_fragments))
    }
    rownames(out) <- names(target_GRanges)
  }

  count_frags <- function(frags) {
    ol <- GenomicRanges::countOverlaps(target_GRanges,
                                       query_fragments[[frags]])
    frag_name <- names(query_fragments)[frags]

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
  }

  if(n_threads > 1) {
    fragment_counts <- mclapply(1:length(query_fragments),
                                count_frags,
                                mc.cores = n_threads)
  } else {
    fragment_counts <- lapply(1:length(query_fragments),
                              count_frags)
  }


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
    out <- out[, names(query_fragments)]
  }

  return(out)
}

#' Filter a GenomicRanges object against another GenomicRanges object
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param target_gr The GenomicRanges object to use for filtering.
#' @param mode Whether to "remove" or "keep" overlaps. Default is "remove"
#' @param ignore_strand Logical, whether or not to ignore the strand of regions in the comparison
#'
#' @return A filtered GenomicRanges object.
#' @export
filter_gr <- function(query_gr,
                      target_gr,
                      mode = "remove",
                      ignore.strand = TRUE) {

  overlapping_fragments <- unique(
    GenomicRanges::queryHits(
      GenomicRanges::findOverlaps(query_gr,
                                  target_gr,
                                  ignore.strand = ignore_strand)))

  if(mode == "remove") {
    query_gr[-overlapping_fragments]
  } else if(mode == "keep") {
    query_gr[overlapping_fragments]
  }
}

#' Filter a list of GenomicRanges objects against a single GenomicRanges object
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param target_gr The GenomicRanges object to use for filtering.
#' @param mode Whether to "remove" or "keep" overlaps. Default is "remove"
#' @param ignore_strand Logical, whether or not to ignore the strand of regions in the comparison
#' @param n_threads A numeric object specifying the number of threads to use. Default is 1.
#'
#' @return A list object containing filtered GenomicRanges objects.
#' @export
filter_fragments <- function(fragment_list,
                             target_gr,
                             mode = "remove",
                             ignore_strand = TRUE,
                             n_threads = 1) {

  if(n_threads == 1) {
    out_list <- lapply(fragment_list,
                       filter_gr,
                       target_gr = target_gr,
                       mode = mode,
                       ignore_strand = ignore_strand)
  } else {
    out_list <- mclapply(fragment_list,
                         filter_gr,
                         target_gr = target_gr,
                         mode = mode,
                         ignore_strand = ignore_strand,
                         mc.cores = n_threads)
  }

  names(out_list) <- names(fragment_list)
  out_list
}
