#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

bam    <- args[1]
bai    <- args[2]
exons  <- args[3]
key    <- args[4]
fasta  <- args[5]
output <- args[6]

# load exons
exons <- readr::read_tsv(exons)
exons <- as.data.frame(exons)

# calculate depth
depth <- ExomeDepth::getBamCounts(
  bed.frame = exons,
  bam = bam,
  index.files = bai,
  referenceFasta = fasta
)

names(depth) <- c('seqnames', 'start', 'end', 'GC', 'counts')
depth$key <- key
depth$exon_id <- exons$exon_id

readr::write_tsv(
  depth,
  output
)
