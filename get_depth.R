#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

bam     <- args[1]
coords  <- args[2]
feature <- args[3]
output  <- args[4]

# load coords
coords <- rtracklayer::import.bed(coords)

# calculate depth
if (feature == "exon") {
  depth <- bamsignals::bamCount(bam, coords)
  gr <- coords
  gr$score <- depth
} else if (feature == "gene") {
  depth <- bamsignals::bamCoverage(bam, coords)
   
  
  gr <- GenomicRanges::tile(coords, width = 1)

  gr <- unlist(gr)
  gr$name <- rep(coords$name, lengths(depth@signals))
  gr$score <- unlist(depth@signals)
} else {
  stop("Invalid feature type. Must be 'exon' or 'gene'.")
}

# write file
rtracklayer::export.bed(
  gr,
  output
)
