context("input/output")

library(ATAComb)


test_that(
  "read_10x_fragments() reads fragments.tsv and generates a list of data.frames, one per cell",
  {
    tenx_dir <- system.file("testdata/outs", package = "ATAComb")

    test_fragments <- read_10x_fragments(outs_dir = tenx_dir,
                                         min_reads = 1000,
                                         remove_chrM = TRUE,
                                         verbose = FALSE)

    expect_true(class(test_fragments) == "list")
    expect_true(length(test_fragments) > 1)
    expect_true("data.frame" %in% class(test_fragments[[1]]))

    fragment_lengths <- unlist(lapply(test_fragments, nrow))

    expect_equal(sum(fragment_lengths > 1000), length(fragment_lengths))

    chrM_fragments <- unlist(lapply(test_fragments, function(x) sum(x$chr == "chrM")))

    expect_equal(sum(chrM_fragments), 0)

  }
)

tenx_dir <- system.file("testdata/outs", package = "ATAComb")

test_fragments <- read_10x_fragments(outs_dir = tenx_dir,
                                     min_reads = 1000,
                                     remove_chrM = TRUE,
                                     verbose = FALSE)

test_that(
  "convert_fragments_gr() converts a list of fragment data.frames to a list of GenomicRanges objects",
  {
    test_gr <- convert_fragments_gr(test_fragments)

    expect_true(class(test_gr) == "list")
    expect_equal(length(test_gr), length(test_fragments))

    test_classes <- unlist(lapply(test_gr, class))

    expect_equal(sum(test_classes == "GRanges"), length(test_classes))

    expect_identical(names(test_gr), names(test_fragments))

    fragment_lengths <- unlist(lapply(test_fragments, nrow))
    gr_lengths <- unlist(lapply(test_gr, length))

    expect_identical(gr_lengths, fragment_lengths)

  }
)


test_that(
  "read_chrom_sizes() reads stored chrom.sizes files and sets n_windows and offsets per chromosome",
  {
    hg38_5e3 <- read_chrom_sizes(genome = "hg38",
                                 window_size = 5e3)

    expect_true("data.frame" %in% class(hg38_5e3))

    window_maxes <- hg38_5e3$n_windows * 5e3
    expect_equal(sum(window_maxes >= hg38_5e3$size), nrow(hg38_5e3))
    expect_equal(sum(window_maxes - 5e3 < hg38_5e3$size), nrow(hg38_5e3))
    expect_equal(sum(diff(hg38_5e3$offset) == hg38_5e3$n_windows[-nrow(hg38_5e3)]), nrow(hg38_5e3) - 1)

    hg19_1e4 <- read_chrom_sizes(genome = "hg19",
                                 window_size = 1e4)

    expect_true("data.frame" %in% class(hg19_1e4))

    window_maxes <- hg19_1e4$n_windows * 1e4
    expect_equal(sum(window_maxes >= hg19_1e4$size), nrow(hg19_1e4))
    expect_equal(sum(window_maxes - 1e4 < hg19_1e4$size), nrow(hg19_1e4))
    expect_equal(sum(diff(hg19_1e4$offset) == hg19_1e4$n_windows[-nrow(hg19_1e4)]), nrow(hg19_1e4) - 1)

    mm10_5e3 <- read_chrom_sizes(genome = "mm10",
                                 window_size = 5e3)

    expect_true("data.frame" %in% class(mm10_5e3))

    window_maxes <- mm10_5e3$n_windows * 5e3
    expect_equal(sum(window_maxes >= mm10_5e3$size), nrow(mm10_5e3))
    expect_equal(sum(window_maxes - 5e3 < mm10_5e3$size), nrow(mm10_5e3))
    expect_equal(sum(diff(mm10_5e3$offset) == mm10_5e3$n_windows[-nrow(mm10_5e3)]), nrow(mm10_5e3) - 1)

    mm9_50 <- read_chrom_sizes(genome = "mm9",
                               window_size = 50)

    expect_true("data.frame" %in% class(mm9_50))

    window_maxes <- mm9_50$n_windows * 50
    expect_equal(sum(window_maxes >= mm9_50$size), nrow(mm9_50))
    expect_equal(sum(window_maxes - 50 < mm9_50$size), nrow(mm9_50))
    expect_equal(sum(diff(mm9_50$offset) == mm9_50$n_windows[-nrow(mm9_50)]), nrow(mm9_50) - 1)
  }
)

hg38_5e3 <- read_chrom_sizes(genome = "hg38",
                             window_size = 5e3)

test_that(
  "fragments_to_window_counts() converts a data.frame of fragments to a table of window counts",
  {

    fragment_df <- test_fragments[[1]]

    frag_windows <- fragments_to_window_counts(fragment_df,
                                               chrom_sizes = hg38_5e3,
                                               window_size = 5e3)

    frag_centers <- (fragment_df$start + fragment_df$end) / 2
    expect_true(frag_centers[1] > as.numeric(names(frag_windows)[1]) * 5e3)
    expect_true(frag_centers[1] < as.numeric(names(frag_windows)[1]) * 5e3 + 5e3)

  }
)

test_that(
  "window_index_to_bed() converts a set of window indexes back to genomic coordinates",
  {

    fragment_df <- test_fragments[[1]]

    frag_windows <- fragments_to_window_counts(fragment_df,
                                               chrom_sizes = hg38_5e3,
                                               window_size = 5e3)

    frag_window_bed <- window_index_to_bed(as.numeric(names(frag_windows)),
                                           chrom_sizes = hg38_5e3,
                                           window_size = 5e3)


  }
)

test_that(
  "convert_fragments_windows() converts a list of fragment data.frames to a sparse matrix of genomic window counts",
  {


  }
)



