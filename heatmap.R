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
col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'exon')
cnv <- cnvr::read_cnv(cnv, col_names)

samples <- stringr::str_split(cnv$sample, '\\.', simplify = TRUE)[, 2]

file_name <- paste(cohort, feature, type, 'png', sep = '.')

png(filename = file_name)

cnvr::plot_heatmap(cnv, column_labels = sort(unique(samples)))

dev.off()
