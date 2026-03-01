#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort    <- args[1]
gene      <- args[2]
key       <- args[3]
obj       <- args[4]

# load cnv object
obj <- readr::read_rds(obj)

# plot
png(paste(cohort, gene, key, 'cnv', 'png', sep = '.'))

exomeCopy::plot.ExomeCopy(
  obj,
  main = paste0(gene, ' (', key, ')')
)

dev.off()
