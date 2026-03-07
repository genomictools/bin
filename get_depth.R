#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

bam         <- args[1]
gene_coords <- args[2]
exon_coords <- args[3]
output      <- args[4]

# count
exon_coords <- rtracklayer::import.bed(exon_coords)
depth <- bamsignals::bamCount(bam, exon_coords)
exon_coords$score <- depth

# write file
rtracklayer::export.bed(
  exon_coords,
  paste(output, 'count', 'bed', sep = ".")
)

# coverage
gene_coords <- rtracklayer::import.bed(gene_coords)
depth <- bamsignals::bamCoverage(bam, gene_coords)

gr <- GenomicRanges::tile(gene_coords, width = 1)
gr <- unlist(gr)
gr$name <- rep(gene_coords$name, lengths(depth@signals))
gr$score <- unlist(depth@signals)

# write file
rtracklayer::export.bed(
  gr,
  paste(output, 'coverage', 'bed', sep = ".")
)
