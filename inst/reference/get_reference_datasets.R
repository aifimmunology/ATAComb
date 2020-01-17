library(data.table)
library(GenomicRanges)

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


