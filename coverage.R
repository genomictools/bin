#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Plot Sequencing Coverage Across Genomic Positions

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args         <- commandArgs(trailingOnly = TRUE)
coverage_bed <- args[1L]
chrom        <- as.character(args[length(args) - 1L])
out_png      <- args[length(args)]

chr_key <- if (startsWith(chrom, "chr")) chrom else paste0("chr", chrom)

cov_gr <- tryCatch(rtracklayer::import.bed(coverage_bed), error = function(e) GenomicRanges::GRanges())
if (length(cov_gr) > 0L) {
  cov_gr <- cov_gr[as.character(GenomicRanges::seqnames(cov_gr)) == chr_key]
}

plot_coverage_tracks(
  cov_gr  = cov_gr,
  chrom   = chr_key,
  out_png = out_png
)

