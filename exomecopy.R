#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(S4Vectors)
})

if (!exists('Rle', mode = 'function')) {
  stop('S4Vectors::Rle is not available in this runtime environment.')
}

cohort      <- args[1]
gene        <- args[2]
test_key    <- args[3]
test_counts <- args[4]
ref_key     <- args[5]
ref_counts  <- args[6]
prob        <- args[7]

# load references
ref_counts <- purrr::map(unlist(strsplit(ref_counts, ',')), rtracklayer::import.bed)
ref_counts <- purrr::map(ref_counts, ~.x$score)
ref_counts <- dplyr::bind_cols(ref_counts)
names(ref_counts) <- unlist(strsplit(ref_key, ','))

# load counts
counts <- rtracklayer::import.bed(test_counts)
names(GenomicRanges::mcols(counts)) <- c('name', test_key)
GenomicRanges::mcols(counts) <- cbind(
  GenomicRanges::mcols(counts),
  ref_counts
)

# generate background
counts$bg <- as.numeric(exomeCopy::generateBackground(unlist(strsplit(ref_key, ',')), counts, median))
counts$log.bg <- log(counts$bg + .1)
counts$width <- GenomicRanges::width(counts)

# call CNVs
library(S4Vectors)
cnv <- exomeCopy::exomeCopy(
  counts,
  test_key,
  # X.names = c("log.bg", "GC", "GC.sq","width"),
  X.names = c("log.bg", "width"),
  goto.cnv = as.numeric(prob)
)

# save cnv object
readr::write_rds(
  cnv,
  paste(cohort, gene, test_key, 'cnv', 'rds', sep = '.')
)

# write output to text file
d <- tibble::as_tibble(exomeCopy::copyCountSegments(cnv))
d <- dplyr::mutate(d, cohort = cohort, gene = gene)
d <- dplyr::rename(d, key = sample.name)

readr::write_tsv(
  d,
  paste(cohort, gene, test_key, 'cnv', 'tsv', sep = '.')
)

# plot
png(paste(cohort, gene, test_key, 'cnv', 'png', sep = '.'))

exomeCopy::plot.ExomeCopy(
  cnv,
  main = paste0(gene, ' (', test_key, ')')
)

dev.off()
