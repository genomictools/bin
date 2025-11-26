#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

signal      <- args[1]
feature     <- args[2]
feature_list<- args[3]
region      <- args[4]
pfb         <- args[5]
output      <- args[6]

# load pfb
pfb <- cnvr::read_pfb(pfb, c('name', 'chr', 'pos', 'pfb'))

# load dbs
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
org <- org.Hs.eg.db::org.Hs.eg.db

# get gene model
if (feature == 'refgene') {
  gene <- unlist(strsplit(feature_list, ','))
  if (length(gene) > 0) {
    gene_models <- cnvr::get_genemodel(txdb, org, gene)
    plot_gene <- TRUE
  } else {
    gene_models <- NULL
    plot_gene <- FALSE
  }
}

# load signal
signal <- cnvr::read_signal(
  signal,
  col_names = c('name', 'baf', 'lrr'),
  pfb = pfb
)

# get overlap
region <- unlist(strsplit(region, '\\_'))
region <- GenomicRanges::GRanges(
  seqnames = region[1],
  ranges = IRanges::IRanges(
    start = as.numeric(region[2]),
    end = as.numeric(region[3])
  )
)

ol <- cnvr::get_overlap(
  region,
  signal,
  flank = GenomicRanges::width(region)/2
)

# LRR
file_name <- paste(output, 'lrr', 'png', sep = '.')
png(filename = file_name, width = 5, height = 5, units = 'in', res = 300)
ylim <- c(c(min(min(ol$lrr), -1)), c(max(max(ol$lrr), 1)))
ylim <- ifelse(ylim > 2, 2, ylim)
ylim <- ifelse(ylim < -2, -2, ylim)

cnvr::plot_signal(
  ol,
  type = 'LRR', ylab = 'LRR',
  ylim = ylim,
  plot_gene = plot_gene,
  gene_model = gene_models
)
dev.off()

# BAF
file_name <- paste(output, 'baf', 'png', sep = '.')
png(filename = file_name, width = 5, height = 5, units = 'in', res = 300)
cnvr::plot_signal(
  ol,
  type = 'BAF', ylab = 'BAF',
  plot_gene = plot_gene,
  gene_model = gene_models
)
dev.off()
