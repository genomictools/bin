#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Ingest Tabular Linkage Output & Generate Visualizations

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args        <- commandArgs(trailingOnly = TRUE)
chrom       <- as.character(args[1L])
in_bed      <- args[2L]
out_png     <- args[3L]
segment_bed <- if (length(args) >= 4L) args[4L] else NULL

lod_gr <- GenomicRanges::GRanges()
if (file.exists(in_bed) && file.info(in_bed)$size > 0L) {
  df <- tryCatch(
    utils::read.table(
      in_bed, header = FALSE, skip = 1L, sep = "\t", stringsAsFactors = FALSE,
      col.names = c("chrom", "start", "end", "name", "lod_single", "lod_multi")
    ),
    error = function(e) NULL
  )
  if (!is.null(df) && nrow(df) > 0L) {
    lod_gr <- GenomicRanges::GRanges(
      seqnames   = as.character(df$chrom),
      ranges     = IRanges::IRanges(start = as.integer(df$start) + 1L, end = as.integer(df$end)),
      name       = df$name,
      lod_single = df$lod_single,
      lod_multi  = df$lod_multi
    )
  }
}

cand_int <- NULL
if (!is.null(segment_bed) && file.exists(segment_bed) && file.info(segment_bed)$size > 0L) {
  seg_gr <- tryCatch(rtracklayer::import.bed(segment_bed), error = function(e) GenomicRanges::GRanges())
  if (length(seg_gr) > 0L) {
    cand_int <- c(GenomicRanges::start(seg_gr[1L]), GenomicRanges::end(seg_gr[1L]))
  }
}

plot_linkage_curve(
  lod_gr             = lod_gr,
  candidate_interval = cand_int,
  chrom              = chrom,
  out_png            = out_png
)
