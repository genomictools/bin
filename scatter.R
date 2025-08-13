#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cnv     <- args[1]
signal  <- args[2]
pfb     <- args[3]
type    <- args[4]
flank   <- args[5]
top_n   <- args[6]

# cnv   <- 'call-cnv-arrays/tests/results/scanned/pheno2.gene.cnv'
# signal<- 'pheno1.sample1.data.txt.adjusted,pheno1.sample3.data.txt.adjusted,pheno2.sample5.data.txt.adjusted,pheno1.sample2.data.txt.adjusted,pheno2.sample4.data.txt.adjusted'
# pfb   <- 'call-cnv-arrays/tests/results/ref/common_all.pfb'
# type <- 'LRR'
# flank <- 1000
# top_n <- 1

# load cnv
d <- readr::read_table(cnv, col_names = FALSE)
d <- dplyr::select(d, region = X1, gene = X9, file = X5)
d <- dplyr::filter(d, gene != 'NOT_FOUND')
d <- transform(d, gene = strsplit(gene, ','))
d <- tidyr::unnest(d, gene)
d <- tidyr::separate(d, file, c('cohort', 'key'), remove = FALSE)
d <- tidyr::separate(d, region, c('chrom', 'start', 'end'), remove = FALSE)
d <- dplyr::group_by(d, gene)
d <- dplyr::reframe(
  d,
  n_variants = length(unique(region)),
  n_samples = length(unique(key)),
  chrom = unique(chrom),
  start = min(as.integer(start)),
  end = max(as.integer(end)),
  cohort = list(unique(cohort)),
  key = list(unique(key)),
  file = list(unique(file))
)

d <- dplyr::top_n(d, top_n, n_samples)

grl <- GenomicRanges::makeGRangesListFromDataFrame(d, split.field = 'gene', keep.extra.columns = TRUE)
# grl <- purrr::map(grl, GenomicRanges::promoters)

# load marker info
pfb <- readr::read_tsv(
  pfb,
  skip = 1,
  col_select = 1:3,
  col_names = c('Name', 'chrom', 'pos')
)

# load signal
signal <- unlist(strsplit(signal, ','))
names(signal) <- signal

signal <- purrr::map_df(signal, readr::read_tsv, .id = 'file')
signal <- dplyr::inner_join(signal, pfb)
signal <- dplyr::mutate(signal, chrom = paste0('chr', chrom))
signal <- GenomicRanges::makeGRangesFromDataFrame(
  signal,
  start.field = 'pos',
  end.field = 'pos', 
  keep.extra.columns = TRUE
)
signal <- GenomicRanges::sort(signal)

# plot
purrr::imap(grl, function(x, .y) {
  gene_name <- .y
  
  sub_signal <- signal[signal$file %in% unlist(x$file)]
  
  ol <- IRanges::subsetByOverlaps(sub_signal, x)
  oll <- split(ol, ol$file)

  purrr::imap(oll, function(x, .y) {
    cohort <- stringr::str_split(.y, '\\.', simplify = TRUE)[1]
    sample <- stringr::str_split(.y, '\\.', simplify = TRUE)[2]
    
    file_name <- paste(cohort, sample, gene_name, type, 'png', sep = '.')
    
    png(filename = file_name, width = 4, height = 4, units = 'in', res = 300)
    
    if ( type == 'lrr' ) {
      plot(
        IRanges::start(x),
        x$`Log R Ratio`,
        ylim = c(-2, 2),
        pch = 19,
        cex = .7,
        xlab = unique(GenomicRanges::seqnames(x)),
        ylab = 'LRR',
        main = paste0(gene_name, " (", cohort, ".", sample, ")")
      )
    } else if ( type == 'baf' ) {
      plot(
        IRanges::start(x),
        x$`B Allele Freq`,
        ylim = c(0, 1),
        pch = 19,
        cex = .7,
        xlab = unique(GenomicRanges::seqnames(x)),
        ylab = 'BAF',
        main = paste0(gene_name, " (", cohort, ".", sample, ")")
      )
    }
    dev.off()
  
  })
})
