#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

gene    <- args[1]
n       <- args[2]
species <- args[3]
genome  <- args[4]
style   <- args[5]
regions <- args[6]

# load gene ids
org <- org.Hs.eg.db::org.Hs.eg.db
if ( gene == 'reference' ) {
  set.seed(1234)
  keys <- AnnotationDbi::keys(org, keytype = 'UNIPROT')
  keys <- sample(keys, as.integer(n))
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

# get exon coordinates for all genes
if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- GenomeInfoDb::keepStandardChromosomes(txdb, species = species, pruning.mode="coarse")

if ( regions == 'null') {
  exons <- GenomicFeatures::exons(txdb, filter = list(gene_id = ids))
} else {
  gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = ids))
  exons <- rtracklayer::import.bed(regions, which = gr)
  GenomicRanges::mcols(exons) <- data.frame(exon_id = exons$name)
}

GenomeInfoDb::seqlevelsStyle(exons) <- style
exons <- as.data.frame(exons)
exons <- dplyr::filter(exons, width > 1)

output <- paste(species, genome, style, gene, 'exons', 'bed', sep = ".")

# write file
readr::write_tsv(
  exons,
  output
)

