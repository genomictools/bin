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
gc <- rtracklayer::import.bed(coords)
gc$score <- CODEX2::getgc(gc, bsgenome)
GenomeInfoDb::seqlevelsStyle(gc) <- style

rtracklayer::export.bed(gc, paste(gene, 'gc.bed', sep = '.'))

# get pmap
pmap <- rtracklayer::import.bed(coords)
pmap$score <- CODEX2::getmapp(pmap, bsgenome)
GenomeInfoDb::seqlevelsStyle(pmap) <- style

rtracklayer::export.bed(pmap, paste(gene, 'pmap.bed', sep = '.'))
