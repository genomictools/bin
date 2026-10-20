#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Disjoin pairwise IBD segments into minimal breakpoint windows

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)

args     <- commandArgs(trailingOnly = TRUE)
ibd_file <- args[1L]
bed_file <- args[2L]

if (!file.exists(ibd_file) || file.info(ibd_file)$size == 0L) {
  file.create(bed_file)
  quit(save = "no", status = 0L)
}

ibd_df <- tryCatch(
  utils::read.table(
    ibd_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("ind1", "hap1", "ind2", "hap2", "chr", "start", "end", "cm")
  ),
  error = function(e) NULL
)

if (is.null(ibd_df) || nrow(ibd_df) == 0L) {
  file.create(bed_file)
  quit(save = "no", status = 0L)
}

ibd_gr <- GenomicRanges::GRanges(
  seqnames = as.character(ibd_df$chr),
  ranges   = IRanges::IRanges(start = as.integer(ibd_df$start), end = as.integer(ibd_df$end)),
  ind1     = as.character(ibd_df$ind1),
  hap1     = as.integer(ibd_df$hap1),
  ind2     = as.character(ibd_df$ind2),
  hap2     = as.integer(ibd_df$hap2),
  cm       = as.numeric(ibd_df$cm)
)

win_gr <- disjoin_segments(ibd_gr)

if (length(win_gr) == 0L) {
  file.create(bed_file)
} else {
  out_df <- data.frame(
    chrom   = as.character(GenomicRanges::seqnames(win_gr)),
    start   = GenomicRanges::start(win_gr),
    end     = GenomicRanges::end(win_gr),
    size_cm = S4Vectors::mcols(win_gr)$size_cm,
    stringsAsFactors = FALSE
  )
  utils::write.table(out_df, file = bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
