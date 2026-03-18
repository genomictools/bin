#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

gene    <- args[1]
genelist<- args[2]
bins    <- args[3]
species <- args[4]
genome  <- args[5]
style   <- args[6]

# load databases
if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- GenomeInfoDb::keepStandardChromosomes(txdb, species = species, pruning.mode="coarse")
org <- org.Hs.eg.db::org.Hs.eg.db

# load gene ids
genelist <- readr::read_lines(genelist)
if (length(genelist) == 0) {
  set.seed(1234)
  keys <- sample(
    AnnotationDbi::keys(org, keytype = 'SYMBOL'),
    as.integer(bins)
  )
} else {
  keys <- genelist
}

ids <- AnnotationDbi::select(
  org,
  keys = keys,
  columns = 'ENTREZID',
  keytype = 'SYMBOL'
)

# get coordinates
for ( feature in c('gene', 'exon')) {
  # get coordinates
  gr <- switch(
    feature,
    'gene' = GenomicFeatures::genes(txdb, filter = list(gene_id = ids$ENTREZID)),
    'exon' = GenomicFeatures::exons(txdb, filter = list(gene_id = ids$ENTREZID), columns = c('gene_id', 'exon_id'))
  )

  # modify object style, remove ranges with width 0
  GenomeInfoDb::seqlevelsStyle(gr) <- style
  gr <- gr[GenomicRanges::width(gr) > 1]

  # edit the name column 
  d <- GenomicRanges::mcols(gr)
  d <- dplyr::mutate_all(as.data.frame(d), ~as.character(purrr::map(.x, ~.x[1])))
  d <- dplyr::left_join(d, ids, by = c('gene_id' = 'ENTREZID'))
  d <- tidyr::unite(d, name, dplyr::everything(), sep = '-')
  GenomicRanges::mcols(gr) <- d

  # write file
  rtracklayer::export.bed(
    gr,
    paste(gene, feature, 'bed', sep = ".")
  )
}
