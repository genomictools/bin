#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Visualize Rare-Variant Drop-Off Scanner & SOR

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args       <- commandArgs(trailingOnly = TRUE)
counts_bed <- args[1L]
seg_bed    <- args[2L]
chrom      <- as.character(args[3L])
png_file   <- args[4L]
cohort     <- if (length(args) >= 5L) args[5L] else ""
family     <- if (length(args) >= 6L) args[6L] else ""

chr_key <- if (startsWith(chrom, "chr")) chrom else paste0("chr", chrom)
chr_num <- sub("^chr", "", chr_key)

sharing_gr <- GenomicRanges::GRanges()
if (file.exists(counts_bed) && file.info(counts_bed)$size > 0L) {
  sh_df <- tryCatch(
    utils::read.table(
      counts_bed, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE,
      col.names = c("chrom", "start", "end", "variant_id", "sharing")
    ),
    error = function(e) NULL
  )
  if (!is.null(sh_df) && nrow(sh_df) > 0L) {
    sh_df <- sh_df[sh_df$chrom != "chrom" & !is.na(as.integer(sh_df$end)), , drop = FALSE]
    keep_chr <- sh_df$chrom == chrom | sh_df$chrom == chr_key | sub("^chr", "", sh_df$chrom) == chr_num
    sh_df <- sh_df[keep_chr, , drop = FALSE]
    if (nrow(sh_df) > 0L) {
      sharing_gr <- GenomicRanges::GRanges(
        seqnames   = rep(chr_key, nrow(sh_df)),
        ranges     = IRanges::IRanges(start = as.integer(sh_df$end), end = as.integer(sh_df$end)),
        variant_id = as.character(sh_df$variant_id),
        sharing    = as.integer(sh_df$sharing)
      )
    }
  }
}

prio_gr <- GenomicRanges::GRanges()
if (file.exists(seg_bed) && file.info(seg_bed)$size > 0L) {
  prio_gr <- tryCatch(rtracklayer::import.bed(seg_bed), error = function(e) GenomicRanges::GRanges())
}

sh_ref <- identify_dropoff(sharing_gr, prio_gr, chrom = chr_key)

plot_dropoff_curve(
  sharing_gr = sh_ref,
  chrom      = chr_key,
  out_png    = png_file,
  cohort     = cohort,
  family     = family
)
