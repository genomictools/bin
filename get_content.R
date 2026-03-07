#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

gene    <- args[1]
coords  <- args[2]
genome  <- args[3]
style   <- args[4]

# select genome
if ( genome == 'hg38') {
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
}

# get gc
library(CODEX2)
gc <- unlist(strsplit(coords, ","))
gc <- purrr::map(gc, rtracklayer::import.bed)
gc <- do.call(c, gc)
gc$score <- CODEX2::getgc(gc, bsgenome)
GenomeInfoDb::seqlevelsStyle(gc) <- style

rtracklayer::export.bed(gc, paste(gene, 'gc.bed', sep = '.'))

# get pmap
pmap <- unlist(strsplit(coords, ","))
pmap <- purrr::map(pmap, rtracklayer::import.bed)
pmap <- do.call(c, pmap)

pmap$score <- CODEX2::getmapp(pmap, bsgenome)
GenomeInfoDb::seqlevelsStyle(pmap) <- style

rtracklayer::export.bed(pmap, paste(gene, 'pmap.bed', sep = '.'))
