library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)

## UCSC Genome Browser chrom.sizes files:
genomes <- c("hg19","hg38","mm10","mm9")

for(genome in genomes) {
  out_file <- paste0("inst/reference/", genome, ".chrom.sizes")
  download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".chrom.sizes"),
                out_file)
}

## UCSC Genome Browser chain files:
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
              "inst/reference/hg38ToHg19.over.chain.gz")
download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
              "inst/reference/hg19ToHg38.over.chain.gz")

## Altius DHS Index from ENCODE:
temp_file <- tempfile(fileext = ".tsv")
download.file("https://www.encodeproject.org/files/ENCFF503GCK/@@download/ENCFF503GCK.tsv",
              temp_file)

alt_idx_raw <- fread(temp_file)

alt_idx_bed <- alt_idx_raw[,c("seqname","start","end","identifier")]
fwrite(alt_idx_bed,
       "inst/reference/hg38_alt_dhs_index.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

alt_idx_gr <- GRanges(seqnames = alt_idx_raw$seqname,
                      IRanges(start = alt_idx_raw$start,
                              end = alt_idx_raw$end),
                      identifier = alt_idx_raw$identifier)

saveRDS(alt_idx_gr,
        "inst/reference/hg38_alt_dhs_index_gr.rds")

## Buenrostro/GSE123577 peak set
temp_file <- tempfile(fileext = ".bed.gz")

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123577&format=file&file=GSE123577%5Fpbmc%5Fpeaks%2Ebed%2Egz",
              temp_file)

temp_chain <- tempfile(fileext = ".over.chain")
R.utils::gunzip("inst/reference/hg19ToHg38.over.chain.gz",
                temp_chain,
                remove = FALSE)

ch <- import.chain(temp_chain)

hg19_peaks <- fread(temp_file)
names(hg19_peaks) <- c("chr","start","end")
hg19_gr <- convert_fragments_gr(list(hg19_peaks))[[1]]

hg38_lo <- liftOver(hg19_gr, ch)
n_lo <- elementNROWS(hg38_lo)

keep_lo <- n_lo == 1

hg38_gr <- unlist(hg38_lo[keep_lo])
hg19_gr <- hg19_gr[keep_lo]

saveRDS(hg19_gr, "inst/reference/hg19_GSE123577_filtered_gr.rds")
saveRDS(hg38_gr, "inst/reference/hg38_GSE123577_filtered_gr.rds")

hg19_bed <- as.data.frame(hg19_gr)[,1:3]
fwrite(hg19_bed,
       "inst/reference/hg19_GSE123577_filtered.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

hg38_bed <- as.data.frame(hg38_gr)[,1:3]
fwrite(hg38_bed,
       "inst/reference/hg38_GSE123577_filtered.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

## TSS regions: ENSEMBL Hg38 v93
## Couldn't download directly, so downloaded through browser from:
## ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/

gtf <- fread("C:/Users/lucasg/Downloads/Homo_sapiens.GRCh38.93.gtf.gz")
names(gtf) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
gtf <- gtf[feature == "gene"]

gtf$attribute <- gsub(" ?[a-z|_]+ \"([A-Za-z0-9|_|-]+)\"",
                      "\\1",
                      gtf$attribute)

gtf <- separate(gtf,
                attribute,
                sep = ";",
                into = c("gene_id","gene_version","gene_name","gene_source","gene_biotype"))

temp_file <- tempfile(fileext = ".h5")
download.file("http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_filtered_feature_bc_matrix.h5",
              temp_file, mode = "wb")

tenx_feat <- H5weaver::read_h5_feature_meta(temp_file)
tenx_ensembl <- tenx_feat$id

keep_gtf <- gtf[gene_id %in% tenx_ensembl]

fwrite(keep_gtf,
       "inst/reference/hg38_ensemble93_tenx_genes.tsv.gz")

gene_gr <- GRanges(seqnames = paste0("chr",gtf$seqname),
                   ranges = IRanges(start = gtf$start,
                           end = gtf$end),
                   strand = gtf$strand,
                   gene_id = gtf$gene_id,
                   gene_name = gtf$gene_name)
saveRDS(gene_gr,
        "inst/reference/hg38_ensemble93_gene_bodies_gr.rds")

gene_bed <- as.data.frame(gene_gr)[,1:4]
fwrite(gene_bed,
       "inst/reference/hg38_ensemble93_gene_bodies.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)

tss_2kb_gr <- resize(gene_gr,
                     width = 4e3,
                     fix = "start")
saveRDS(tss_2kb_gr,
        "inst/reference/hg38_ensemble93_tss_2kb_gr.rds")

tss_2kb_bed <- as.data.frame(tss_2kb_gr)[,1:4]
fwrite(tss_2kb_bed,
       "inst/reference/hg38_ensemble93_tss_2kb.bed.gz",
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)
