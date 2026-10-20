#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Prioritize Shared Haplotype Segments Across a Family Pedigree

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)
library(pedtools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 8L) {
  windows_file      <- args[1L]
  ibd_file          <- args[2L]
  ped_file          <- args[3L]
  pheno_file        <- args[4L]
  allow_mismatch    <- as.integer(args[5L])
  allow_phenocopies <- as.integer(args[6L])
  out_file          <- args[7L]
  out_segment_bed   <- args[8L]
} else if (length(args) == 7L) {
  windows_file      <- args[1L]
  ibd_file          <- args[2L]
  ped_file          <- args[3L]
  pheno_file        <- args[4L]
  allow_mismatch    <- as.integer(args[5L])
  allow_phenocopies <- as.integer(args[6L])
  out_file          <- args[7L]
  out_segment_bed   <- sub("\\.bed$", ".segment.bed", out_file)
} else {
  windows_file      <- args[1L]
  ibd_file          <- args[2L]
  ped_file          <- args[3L]
  pheno_file        <- NULL
  allow_mismatch    <- as.integer(args[4L])
  allow_phenocopies <- as.integer(args[5L])
  out_file          <- args[6L]
  out_segment_bed   <- if (length(args) >= 7L) args[7L] else sub("\\.bed$", ".segment.bed", out_file)
}

header_cols <- c(
  "#chrom", "start", "end", "score", "size_mb", "size_cm", "snps",
  "individuals", "aff", "obligate_carriers", "unaff",
  "strand_mismatches", "transmission_breaks", "rank"
)

write_prioritized_bed <- function(gr, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(header_cols, collapse = "\t"), con)
  if (length(gr) > 0L) {
    mc <- S4Vectors::mcols(gr)
    df <- data.frame(
      chrom               = as.character(GenomicRanges::seqnames(gr)),
      start               = GenomicRanges::start(gr) - 1L,
      end                 = GenomicRanges::end(gr),
      score               = mc$score,
      size_mb             = mc$size_mb,
      size_cm             = mc$size_cm,
      snps                = mc$snps,
      individuals         = mc$individuals,
      aff                 = mc$aff,
      obligate_carriers   = mc$obligate_carriers,
      unaff               = mc$unaff,
      strand_mismatches   = mc$strand_mismatches,
      transmission_breaks = mc$transmission_breaks,
      rank                = mc$rank,
      stringsAsFactors    = FALSE
    )
    utils::write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

if (!file.exists(windows_file) || file.info(windows_file)$size == 0L ||
    !file.exists(ibd_file) || file.info(ibd_file)$size == 0L) {
  write_prioritized_bed(GenomicRanges::GRanges(), out_file)
  file.create(out_segment_bed)
  quit(save = "no", status = 0L)
}

# Ingest windows file (supporting both 4/5-col table and 6/7-col BED format)
win_df <- tryCatch(
  utils::read.table(windows_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
  error = function(e) NULL
)
if (is.null(win_df) || nrow(win_df) == 0L) {
  write_prioritized_bed(GenomicRanges::GRanges(), out_file)
  file.create(out_segment_bed)
  quit(save = "no", status = 0L)
}

if (ncol(win_df) >= 7L) {
  # 6-col BED + annot-tsv count
  w_chrom   <- as.character(win_df[[1L]])
  w_start   <- as.integer(win_df[[2L]]) + 1L
  w_end     <- as.integer(win_df[[3L]])
  w_size_cm <- suppressWarnings(as.numeric(win_df[[5L]]))
  if (all(is.na(w_size_cm) | w_size_cm == 0)) w_size_cm <- (w_end - w_start + 1L) / 1e6
  w_snps    <- as.integer(win_df[[7L]])
} else {
  # 4/5-col table: chrom, start, end, size_cm [, snps]
  w_chrom   <- as.character(win_df[[1L]])
  w_start   <- as.integer(win_df[[2L]])
  w_end     <- as.integer(win_df[[3L]])
  w_size_cm <- suppressWarnings(as.numeric(win_df[[4L]]))
  if (all(is.na(w_size_cm) | w_size_cm == 0)) w_size_cm <- (w_end - w_start + 1L) / 1e6
  w_snps    <- if (ncol(win_df) >= 5L) as.integer(win_df[[5L]]) else rep(0L, nrow(win_df))
}

win_gr <- GenomicRanges::GRanges(
  seqnames = w_chrom,
  ranges   = IRanges::IRanges(start = w_start, end = w_end),
  size_cm  = w_size_cm,
  snps     = w_snps
)

# Ingest IBD table
ibd_df <- tryCatch(
  utils::read.table(
    ibd_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("ind1", "hap1", "ind2", "hap2", "chr", "start", "end", "cm")
  ),
  error = function(e) NULL
)
if (is.null(ibd_df) || nrow(ibd_df) == 0L) {
  write_prioritized_bed(GenomicRanges::GRanges(), out_file)
  file.create(out_segment_bed)
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

# Load Pedigree using pedtools
ped <- tryCatch(pedtools::readPed(ped_file), error = function(e) NULL)
if (is.null(ped)) {
  ped_df <- utils::read.table(ped_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ped_df) <- tolower(colnames(ped_df))
  if (ncol(ped_df) >= 5L && !("id" %in% colnames(ped_df))) {
    colnames(ped_df)[1:5] <- c("famid", "id", "fid", "mid", "sex")
  }
  ped <- pedtools::ped(id = as.character(ped_df$id), fid = as.character(ped_df$fid), mid = as.character(ped_df$mid), sex = as.integer(ped_df$sex))
}

status <- if (!is.null(pheno_file) && nzchar(pheno_file) && file.exists(pheno_file) && file.info(pheno_file)$size > 0L) {
  ph_df <- utils::read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ph_df) <- tolower(colnames(ph_df))
  affected <- as.character(ph_df[[1L]][ph_df[[2L]] == 2])
  infer_lineage(ped, affected = affected)
} else if ("aff" %in% colnames(ped_df)) {
  if (any(ped_df$aff == 3)) {
    stats::setNames(as.integer(ped_df$aff), as.character(ped_df$id))
  } else {
    infer_lineage(ped, affected = as.character(ped_df$id[ped_df$aff == 2]))
  }
} else {
  infer_lineage(ped, affected = character(0))
}

prio_gr <- prioritize_segments(
  windows           = win_gr,
  ibd               = ibd_gr,
  ped               = ped,
  status            = status,
  allow_mismatch    = allow_mismatch,
  allow_phenocopies = allow_phenocopies
)

write_prioritized_bed(prio_gr, out_file)

if (length(prio_gr) > 0L) {
  rtracklayer::export.bed(prio_gr[1L], out_segment_bed)
} else {
  file.create(out_segment_bed)
}
