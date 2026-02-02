#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort    <- args[1]
gene      <- args[2]
key       <- args[3]
obj       <- args[4]
threshold <- args[5]

# load cnv object
cnv <- readr::read_rds(obj)

# plot
png(paste(cohort, gene, key, 'cnv', 'png', sep = '.'))
anno <- cnv@annotations
flank <- 10000

ExomeDepth::plot (
  cnv,
  sequence = unique(as.character(anno$chromosome)),
  xlim = c(min(anno$start) - flank, max(anno$end) + flank),
  main = paste0(gene, ' (', key, ')'),
  count.threshold = as.integer(threshold),
  cex.lab = 0.8,
  with.gene = TRUE,
  xlab = unique(as.character(anno$chromosome))
)

dev.off()
