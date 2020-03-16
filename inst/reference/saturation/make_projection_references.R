library(purrr)
library(ATAComb)
options(stringsAsFactors = FALSE)

dataset_names <- c("tenx_nextgem",
                   "leukopak",
                   "high_neutrophil")

summary_csv_files <- list.files("inst/reference/saturation",
                                pattern = "_summary.csv",
                                full.names = TRUE)
saturation_files <- list.files("inst/reference/saturation",
                               pattern = "_saturation_curve.tsv.gz",
                               full.names = TRUE)

summaries <- map(summary_csv_files, read.csv)
saturations <- map(saturation_files, read.table, header = F, col.names = c("count", "freq"))

saturations <- map(saturations, as.matrix)

metrics <- map(1:length(summaries),
               function(x) {
                 summary <- summaries[[x]]
                 saturation <- saturations[[x]]
                 data.frame(total_reads = summary$num_fragments,
                            total_umis = summary$total_usable_fragments,
                            total_counts = sum(saturation[,"count"] * saturation[,"freq"]))
               })

projections <- map(1:length(summaries),
                   function(x) {
                     saturation <- saturations[[x]]
                     metric <- metrics[[x]]
                     suppressWarnings(diversity_projection(saturation, metric, max_val = 2e9))
                   })

results <- map(1:length(projections),
               function(x) {
                 projection <- projections[[x]]
                 projection$dataset <- dataset_names[x]
                 projection
               })
results <- do.call(rbind, results)

write.csv(results,
          "inst/reference/saturation/reference_projections.csv")

library(ggplot2)

ggplot() +
  geom_line(data = results,
            aes(x = n_raw_reads,
                y = expected_umis,
                group = dataset,
                color = dataset))
