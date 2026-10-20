#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Sum Singlepoint & Multipoint Parametric Linkage LOD Scores

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args            <- commandArgs(trailingOnly = TRUE)
n_args          <- length(args)
markers_file    <- args[1L]
min_linkage     <- as.numeric(args[n_args])
out_segment_bed <- args[n_args - 1L]
out_bed         <- args[n_args - 2L]
chrom           <- as.character(args[n_args - 3L])
tbl_files       <- args[seq_len(n_args - 4L)[-1L]]

lod_gr <- sum_linkage_tables(
  markers_file = markers_file,
  tbl_files    = tbl_files,
  chrom        = chrom,
  min_linkage  = min_linkage
)

if (length(lod_gr) == 0L) {
  file.create(out_bed)
  file.create(out_segment_bed)
  quit(save = "no", status = 0L)
}

# Export detailed linkage LOD score table
con <- file(out_bed, open = "wt")
writeLines("#chrom\tstart\tend\tname\tlod_single\tlod_multi", con)
df <- data.frame(
  chrom      = as.character(GenomicRanges::seqnames(lod_gr)),
  start      = GenomicRanges::start(lod_gr) - 1L,
  end        = GenomicRanges::end(lod_gr),
  name       = S4Vectors::mcols(lod_gr)$name,
  lod_single = S4Vectors::mcols(lod_gr)$lod_single,
  lod_multi  = S4Vectors::mcols(lod_gr)$lod_multi,
  stringsAsFactors = FALSE
)
utils::write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

# Export top recommended candidate interval if found
c_int <- attr(lod_gr, "candidate_interval")
if (!is.null(c_int) && all(is.finite(c_int))) {
  peak_val <- max(S4Vectors::mcols(lod_gr)$lod_multi, na.rm = TRUE)
  top_gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = as.integer(c_int[1L]), end = as.integer(c_int[2L])),
    score    = if (is.finite(peak_val)) round(peak_val, 4L) else 0
  )
  rtracklayer::export.bed(top_gr, out_segment_bed)
} else {
  file.create(out_segment_bed)
}
