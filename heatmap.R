#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort    <- args[1]
tool      <- args[2]
feature   <- args[3]
cnv       <- args[4]
type      <- args[5]
genelist  <- args[6]
n_samples <- args[7]
n_genes   <- args[8]

# load cnv
col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'distance')
cnv <- cnvr::read_cnv(cnv, col_names)

genelist <- unlist(strsplit(genelist, ','))

if (length(genelist) > 1) {
  cnv$gene <- purrr::map_chr(strsplit(cnv$gene, ','), ~paste(intersect(.x, genelist), collapse = ','))
  cnv <- cnv[cnv$gene != '']
}

file_name <- paste(cohort, tool, feature, type, 'heatmap', 'png', sep = '.')

png(filename = file_name)

cnvr::plot_heatmap(
    cnv,
    top_samples = as.integer(n_samples),
    top_genes = as.integer(n_genes)
)

dev.off()
