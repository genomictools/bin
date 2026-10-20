#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Autonomous Rare-Variant Drop-Off Scanner & SOR Identification

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args            <- commandArgs(trailingOnly = TRUE)
counts_bed      <- args[1L]
chrom           <- as.character(args[2L])
out_bed         <- args[3L]
out_segment_bed <- if (length(args) >= 4L) args[4L] else sub("\\.bed$", ".segment.bed", out_bed)
prioritized_bed <- if (length(args) >= 5L && nzchar(args[5L])) args[5L] else NULL

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
if (!is.null(prioritized_bed) && file.exists(prioritized_bed) && file.info(prioritized_bed)$size > 0L) {
  prio_gr <- tryCatch(rtracklayer::import.bed(prioritized_bed), error = function(e) GenomicRanges::GRanges())
}

sh_ref     <- identify_dropoff(sharing_gr, prio_gr, chrom = chr_key)
rgn1       <- attr(sh_ref, "rgn1")
dropoff_df <- attr(sh_ref, "dropoff_df")

if (!is.null(dropoff_df) && nrow(dropoff_df) > 0L) {
  colnames(dropoff_df)[1L] <- "#chrom"
  utils::write.table(dropoff_df, file = out_bed, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
  cat("#chrom\tstart\tend\tregion\tcarrier_count\tn_variants\tsize_mb\tdescription\n", file = out_bed)
}

if (!is.null(rgn1) && !any(is.na(rgn1))) {
  top_gr <- GenomicRanges::GRanges(
    seqnames = chr_key,
    ranges   = IRanges::IRanges(start = as.integer(rgn1[1L]), end = as.integer(rgn1[2L]))
  )
  rtracklayer::export.bed(top_gr, out_segment_bed)
} else {
  file.create(out_segment_bed)
}
