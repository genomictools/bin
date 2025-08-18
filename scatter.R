#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cnv     <- args[1]
signal  <- args[2]
pfb     <- args[3]
type    <- args[4]
flank   <- args[5]
top_n   <- args[6]

# load cnv
col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'exon')
cnv <- cnvr::read_cnv(cnv, col_names)
cnv <- cnv[cnv$gene != 'NOT_FOUND']
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
        # get gene model
        gene <- unlist(strsplit(.x$gene, ','))
        if ( length(gene) ) {
          gene <- intersect(gene, keys)
          gene_models <- cnvr::get_genemodel(txdb, org, gene)

          # get overlap
          ol <- cnvr::get_overlap(
            .x,
            signal,
            flank = GenomicRanges::width(.x)/2
          )

          # make plot
          cohort <- unlist(strsplit(.x$sample, '\\.'))[1]
          sample <- unlist(strsplit(.x$sample, '\\.'))[2]
          file_name <- paste(cohort, sample, .x$region, type, 'png', sep = '.')
          png(filename = file_name, width = 5, height = 4 + length(gene) / 2, units = 'in', res = 300)
          cnvr::plot_signal(
            ol,
            type = toupper(type), ylab = toupper(type),
            gene_model = gene_models
          )
          dev.off()
        }
      }
    )
  }
)
