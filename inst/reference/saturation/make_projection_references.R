library(purrr)
library(ATAComb)
options(stringsAsFactors = FALSE)

summary_csv_files <- list.files("inst/reference/saturation",
                                pattern = "_summary.csv",
                                full.names = TRUE)
saturation_files <- list.files("inst/reference/saturation",
                               pattern = "_saturation_curve.tsv.gz",
                               full.names = TRUE)

summaries <- map(summary_csv_files, read.csv)
saturations <- map(saturation_files, read.table, header = F, col.names = c("count", "freq"))

saturations <- map(saturations, as.matrix)

total_metrics <- data.frame(total_reads = metrics_summary$num_fragments,
                            total_umis = metrics_summary$total_usable_fragments,
                            total_counts = sum(saturation_mat[,"count"] * saturation_mat[,"frequency"]))

saturation_projection <- suppressWarnings(
  diversity_projection(saturation_mat,
                       total_metrics,
                       max_val = 2e9)
)
