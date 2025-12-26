#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

signal      <- args[1]
feature     <- args[2]
feature_list<- args[3]
region      <- args[4]
pfb         <- args[5]
cytobands   <- args[6]
output      <- args[7]

# load pfb
pfb <- cnvr::read_pfb(pfb, c('name', 'chr', 'pos', 'pfb'))

# load region
region <- unlist(strsplit(region, '\\_'))
region <- GenomicRanges::GRanges(
  seqnames = region[1],
  ranges = IRanges::IRanges(
    start = as.numeric(region[2]),
    end = as.numeric(region[3])
  )
)

flank <- GenomicRanges::width(region)

# load dbs
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
org <- org.Hs.eg.db::org.Hs.eg.db

# get gene list
if (feature == 'refgene') {
  gene <- unlist(strsplit(feature_list, ','))
} else {
  gene_id <- unique(IRanges::subsetByOverlaps(GenomicFeatures::genes(txdb), region)$gene_id)
  gene <- AnnotationDbi::select(org, gene_id, 'SYMBOL', 'ENTREZID')$SYMBOL  
}

# get gene model
if (length(gene) > 0) {
  gene_models <- cnvr::get_genemodel(txdb, org, gene)
} else {
  gene_models <- NULL
}

# load signal
signal <- cnvr::read_signal(
  signal,
  col_names = c('name', 'baf', 'lrr'),
  pfb = pfb
)

# read cytobands
col_names <- c('num', 'chrom', 'start', 'end', 'band', 'name', 'stain')
cytobands <- cnvr::read_cytobands(cytobands, col_names = col_names)

# LRR
file_name <- paste(output, 'lrr', 'png', sep = '.')
height <- 4 + length(gene) / 1.5
png(filename = file_name, width = 5, height = height, units = 'in', res = 300)

cnvr::plot_signal(
  signal, region,
  flank = flank,
  type = 'LRR', ylab = 'LRR',
  ylim = c(-1.5, 1.5),
  gene_model = gene_models,
  bands = cytobands
)

dev.off()

# BAF
file_name <- paste(output, 'baf', 'png', sep = '.')
png(filename = file_name, width = 5, height = 5, units = 'in', res = 300)
cnvr::plot_signal(
  flank = flank,
  signal, region,
  type = 'BAF', ylab = 'BAF',
  gene_model = gene_models
)
dev.off()
