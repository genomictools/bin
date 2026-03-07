#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

gene    <- args[1]
# feature <- args[2]
bins    <- args[2]
species <- args[3]
genome  <- args[4]
style   <- args[5]
# output  <- args[6]

# load gene ids
org <- org.Hs.eg.db::org.Hs.eg.db
if ( gene == 'reference' ) {
  set.seed(1234)
  keys <- AnnotationDbi::keys(org, keytype = 'UNIPROT')
  keys <- sample(keys, as.integer(bins))
  ids <- AnnotationDbi::select(
    org,
    keys = keys,
    columns = 'ENTREZID',
    keytype = 'UNIPROT'
  )
} else {
  keys <- unlist(strsplit(gene, split = ","))
  ids <- AnnotationDbi::select(
    org,
    keys = keys,
    columns = 'ENTREZID',
    keytype = 'SYMBOL'
  )
}

ids <- unique(dplyr::pull(ids, ENTREZID))

# get coordinates
if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- GenomeInfoDb::keepStandardChromosomes(txdb, species = species, pruning.mode="coarse")

# if ( feature == 'exon' ) {
#   coords <- GenomicFeatures::exons(txdb, filter = list(gene_id = ids))
# } else if ( feature == 'gene' ) {
#   coords <- GenomicFeatures::genes(txdb, filter = list(gene_id = ids))
# }

# gene
genes <- GenomicFeatures::genes(txdb, filter = list(gene_id = ids))
GenomeInfoDb::seqlevelsStyle(genes) <- style
names(GenomicRanges::mcols(genes)) <- 'name'
genes <- genes[GenomicRanges::width(genes) > 1]

# write file
rtracklayer::export.bed(
  genes,
  paste(gene, "gene", 'bed', sep = ".")
)

# exon
exon <- GenomicFeatures::exons(txdb, filter = list(gene_id = ids))

GenomeInfoDb::seqlevelsStyle(exon) <- style
names(GenomicRanges::mcols(exon)) <- 'name'
exon <- exon[GenomicRanges::width(exon) > 1]

# write file
rtracklayer::export.bed(
  exon,
  paste(gene, "exon", 'bed', sep = ".")
)
