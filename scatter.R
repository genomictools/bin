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
cnv <- split(cnv, cnv$region)

# load pfb
pfb <- cnvr::read_pfb(pfb, c('name', 'chr', 'pos', 'pfb'))

# load signal
col_names <- c('name', 'baf', 'lrr')
signal <- unlist(strsplit(signal, ','))
names(signal) <- signal

signal <- purrr::imap(
  signal,
  ~{
    gr <- cnvr::read_signal(.x, col_names = col_names, pfb = pfb)
    gr$sample <- .y
    gr
  }
)
signal <- GenomicRanges::GRangesList(signal)
signal <- unlist(signal)

# load dbs
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
org <- org.Hs.eg.db::org.Hs.eg.db
keys <- AnnotationDbi::keys(org, 'SYMBOL')

purrr::imap(
  cnv,
  ~{
    # get gene model
    gene <- unlist(strsplit(.x$gene, ','))
    gene <- head(intersect(gene, keys), 2)
    gene_models <- cnvr::get_genemodel(txdb, org, gene)
    
    # get overlap
    sub_signal <- signal[signal$sample == .x$sample]
    ol <- cnvr::get_overlap(
      .x,
      sub_signal,
      flank = GenomicRanges::width(.x)/2
    )
    
    # make plot
    cohort <- unlist(strsplit(.x$sample, '\\.'))[1]
    sample <- unlist(strsplit(.x$sample, '\\.'))[2]
    file_name <- paste(cohort, sample, .x$region, type, 'png', sep = '.')
    png(filename = file_name, width = 4, height = 4, units = 'in', res = 300)
    cnvr::plot_signal(
      ol,
      type = toupper(type), ylab = toupper(type),
      gene_model = gene_models
    )
    dev.off()
  }
)
