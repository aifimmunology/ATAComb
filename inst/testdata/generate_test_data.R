library(ATAComb)
library(data.table)

outs_dir <- "C:/datasets/atac_perm_cells/full_tenx_run/perm_cell_dig_0.01_analysis"

fragments_file <- file.path(outs_dir, "fragments.tsv.gz")
singlecell_file <- file.path(outs_dir, "singlecell.csv")

singlecell <- data.table::fread(singlecell_file)
singlecell <- singlecell[-1,]

min_reads <- 1000

set.seed(3030)
singlecell_pass_cell <- singlecell[sample(which(passed_filters >= min_reads & is__cell_barcode == 1), 50),]
singlecell_pass_noncell <- singlecell[sample(which(passed_filters >= min_reads & is__cell_barcode != 1), 50),]
singlecell_fail <- singlecell[sample(which(passed_filters < min_reads), 100),]

singlecell <- rbind(singlecell_pass_cell,
                    singlecell_pass_noncell,
                    singlecell_fail)

fwrite(singlecell,
       "inst/testdata/outs/singlecell.csv")

fragments <- vroom::vroom(fragments_file,
                          delim = "\t",
                          col_names = c("chr","start","end","barcode","n_reads"),
                          col_types = c(chr = "c", start = "i", end = "i", barcode = "c", n_reads = "i"))

fragments <- fragments[fragments$barcode %in% singlecell$barcode,]

fwrite(fragments,
       "inst/testdata/outs/fragments.tsv.gz",
       col.names = FALSE,
       sep = "\t")

file.copy(file.path(outs_dir,"summary.csv"),
          "inst/testdata/outs/summary.csv")

