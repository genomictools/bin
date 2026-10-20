#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Visualize IBD Haplotypes and Minimal Shared Region (SOR)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)
library(pedtools)

args              <- commandArgs(trailingOnly = TRUE)
ibd_file          <- args[1L]
ped_file          <- args[2L]
pheno_file        <- args[3L]
allow_phenocopies <- if (length(args) >= 4L) as.integer(args[4L]) else 2L
png_file          <- if (length(args) >= 5L) args[5L] else "haplotypes.png"
seg_file          <- if (length(args) >= 6L && nzchar(args[6L])) args[6L] else NULL

# Ingest IBD segments
ibd_df <- tryCatch(
  utils::read.table(
    ibd_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
    col.names = c("ind1", "hap1", "ind2", "hap2", "chr", "start", "end", "cm")
  ),
  error = function(e) NULL
)

if (is.null(ibd_df) || nrow(ibd_df) == 0L) {
  plot_haplotype_matrix(matrix(numeric(0), 0, 0), out_png = png_file)
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

# Load pedigree using pedtools
ped <- tryCatch(pedtools::readPed(ped_file), error = function(e) NULL)
if (is.null(ped)) {
  ped_df <- utils::read.table(ped_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ped_df) <- tolower(colnames(ped_df))
  if (ncol(ped_df) >= 5L && !("id" %in% colnames(ped_df))) {
    colnames(ped_df)[1:5] <- c("famid", "id", "fid", "mid", "sex")
  }
  ped <- pedtools::ped(id = as.character(ped_df$id), fid = as.character(ped_df$fid), mid = as.character(ped_df$mid), sex = as.integer(ped_df$sex))
}

affected <- character(0)
if (!is.null(pheno_file) && nzchar(pheno_file) && file.exists(pheno_file) && file.info(pheno_file)$size > 0L) {
  ph_df <- utils::read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ph_df) <- tolower(colnames(ph_df))
  affected <- as.character(ph_df[[1L]][ph_df[[2L]] == 2])
} else if ("aff" %in% colnames(ped_df)) {
  affected <- as.character(ped_df$id[ped_df$aff == 2])
}

status   <- infer_lineage(ped, affected = affected)
carriers <- names(status)[status %in% c(2L, 3L)]

m <- tile_ibd_matrix(
  ibd               = ibd_gr,
  carriers          = carriers,
  tile_width        = 1000000L,
  allow_phenocopies = allow_phenocopies
)

# If recommended candidate segment is supplied, override shared columns
if (!is.null(seg_file) && file.exists(seg_file) && file.info(seg_file)$size > 0L) {
  seg_df <- tryCatch(
    utils::read.table(seg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(seg_df) && nrow(seg_df) > 0L && ncol(seg_df) >= 3L && ncol(m) > 0L) {
    seg_start  <- as.numeric(seg_df[[2L]][1L])
    seg_end    <- as.numeric(seg_df[[3L]][1L])
    act_tiles  <- as.integer(colnames(m))
    start_tile <- max(1L, floor(seg_start / 1e6) + 1L)
    end_tile   <- ceiling(seg_end / 1e6)
    match_cols <- which(act_tiles >= start_tile & act_tiles <= end_tile)
    if (length(match_cols) > 0L) {
      attr(m, "shared_cols") <- match_cols
    }
  }
}

plot_haplotype_matrix(m, out_png = png_file)
