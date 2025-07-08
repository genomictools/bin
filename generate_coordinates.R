#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

chrom   <- args[1]
start   <- args[2]
end     <- args[3]
genome  <- args[4]
style   <- args[5]
coding  <- args[6]
chunk   <- args[7]
output  <- args[8]

# load genes
if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
if (genome == 'hg38') bsdb <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

if ( start != 'null' && end != 'null') {
  q <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = as.integer(start), end = as.integer(end))
  )
} else {
  q <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(bsdb)[chrom])
  )
}

if ( coding == 'true' ) {
  gene_coordinates <- GenomicFeatures::genes(
    txdb,
    filter = list(tx_chrom = chrom),
    columns = AnnotationDbi::columns(txdb)
  )
  gene_coordinates <- IRanges::subsetByOverlaps(gene_coordinates, q)
} else if ( coding == 'false' ) {
  gene_coordinates <- unlist(GenomicRanges::tile(q, width = as.integer(chunk)))
} else {
  stop("coding can be 'true' or 'false'.")
}

GenomeInfoDb::seqlevels(gene_coordinates) <- chrom
GenomeInfoDb::seqlevelsStyle(gene_coordinates) <- style

rtracklayer::export.bed(
  gene_coordinates,
  output
)
