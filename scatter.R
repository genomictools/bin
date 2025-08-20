#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort  <- args[1]
gene    <- args[2]
cnv     <- args[3]
signal  <- args[4]
pfb     <- args[5]
flank   <- args[6]

# load cnv
col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf')
cnv <- cnvr::read_cnv(cnv, col_names)
cnv <- split(cnv, cnv$sample)

# load pfb
pfb <- cnvr::read_pfb(pfb, c('name', 'chr', 'pos', 'pfb'))

# load dbs
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
org <- org.Hs.eg.db::org.Hs.eg.db

txdb_keys <- AnnotationDbi::keys(txdb, 'GENEID')
org_keys <- AnnotationDbi::keys(org, 'ENTREZID')

keys <- intersect(txdb_keys, org_keys)
keys <- AnnotationDbi::select(org, keys, 'SYMBOL', 'ENTREZID')$SYMBOL

# get gene model
if (length(intersect(gene, keys)) > 0) {
  gene_models <- cnvr::get_genemodel(txdb, org, gene)
  plot_gene <- TRUE
} else {
  gene_models <- NULL
  plot_gene <- FALSE
}

purrr::imap(
  cnv,
  ~{
    # load signal
    signal <- cnvr::read_signal(
      unique(.x$sample),
      col_names = c('name', 'baf', 'lrr'),
      pfb = pfb
    )

    # split by region
    gr <- split(.x, .x$region)
    purrr::imap(
      gr,
      ~{
        # get overlap
        ol <- cnvr::get_overlap(
          .x,
          signal,
          flank = GenomicRanges::width(.x)/2
        )

        # make plot
        cohort <- unlist(strsplit(.x$sample, '\\.'))[1]
        sample <- unlist(strsplit(.x$sample, '\\.'))[2]

        # LRR
        file_name <- paste(cohort, gene, sample, .x$region, 'lrr', 'png', sep = '.')
        png(filename = file_name, width = 4, height = 4, units = 'in', res = 300)
        cnvr::plot_signal(
          ol,
          type = 'LRR', ylab = 'LRR',
          plot_gene = plot_gene,
          gene_model = gene_models
        )
        dev.off()
        
        # BAF
        file_name <- paste(cohort, gene, sample, .x$region, 'baf', 'png', sep = '.')
        png(filename = file_name, width = 4, height = 4, units = 'in', res = 300)
        cnvr::plot_signal(
          ol,
          type = 'BAF', ylab = 'BAF',
          plot_gene = plot_gene,
          gene_model = gene_models
        )
        dev.off()
      }
    )
  }
)
