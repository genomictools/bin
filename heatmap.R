#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort    <- args[1]
feature   <- args[2]
cnv       <- args[3]
type      <- args[4]
n_samples <- args[5]
n_genes   <- args[6]

# load cnv
col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'distance')
cnv <- cnvr::read_cnv(cnv, col_names)

file_name <- paste(cohort, feature, type, 'heatmap', 'png', sep = '.')

png(filename = file_name)

cnvr::plot_heatmap(
    cnv,
    top_samples = as.integer(n_samples),
    top_genes = as.integer(n_genes)
)

dev.off()
