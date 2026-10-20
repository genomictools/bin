#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Prepare pedigree branches and marker files for MERLIN parametric linkage

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(pedtools)

args       <- commandArgs(trailingOnly = TRUE)
vcf_file   <- args[1L]
ped_file   <- args[2L]
pheno_file <- args[3L]
chrom       <- as.character(args[4L])
out_dir     <- args[5L]
max_markers <- if (length(args) >= 6L && nzchar(args[6L]) && args[6L] != "null") as.numeric(args[6L]) else 120L

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

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

# Determine dynamic family ID
famid_val <- tryCatch(as.character(ped$FAMID)[1L], error = function(e) "")
if (!nzchar(famid_val)) {
  famid_val <- tryCatch(as.character(ped_df$famid)[1L], error = function(e) "fam01")
}
if (!nzchar(famid_val)) famid_val <- "fam01"

affected <- character(0)
if (file.exists(pheno_file) && file.info(pheno_file)$size > 0L) {
  ph_df <- utils::read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ph_df) <- tolower(colnames(ph_df))
  affected <- as.character(ph_df[[1L]][ph_df[[2L]] == 2])
} else if (exists("ped_df") && "aff" %in% colnames(ped_df)) {
  affected <- as.character(ped_df$id[ped_df$aff == 2])
}

status   <- infer_lineage(ped, affected = affected)
carriers <- names(status)[status %in% c(2L, 3L)]

vcf <- tryCatch(suppressWarnings(VariantAnnotation::readVcf(vcf_file)), error = function(e) NULL)
if (is.null(vcf) || nrow(vcf) == 0L) {
  cat("chrom\tmarker\tpos\taf\tname\n", file = file.path(out_dir, "markers.tsv"))
  cat("n_branches\tn_markers\n0\t0\n", file = file.path(out_dir, "counts.tsv"))
  quit(save = "no", status = 0L)
}

gt      <- VariantAnnotation::geno(vcf)$GT
rr      <- SummarizedExperiment::rowRanges(vcf)
pos     <- as.integer(GenomicRanges::start(rr))
info_df <- VariantAnnotation::info(vcf)

af_vec <- rep(NA_real_, nrow(gt))
if ("gnomADg_AF_nfe" %in% colnames(info_df)) {
  raw_gn <- as.list(info_df$gnomADg_AF_nfe)
  af_vec <- vapply(raw_gn, function(x) {
    if (length(x) >= 1L) suppressWarnings(as.numeric(x[1L])) else NA_real_
  }, numeric(1))
}
if ("CAF" %in% colnames(info_df) && any(is.na(af_vec))) {
  raw_caf <- as.list(info_df$CAF)
  caf_af  <- vapply(raw_caf, function(x) {
    if (length(x) >= 2L) suppressWarnings(as.numeric(x[2L])) else NA_real_
  }, numeric(1))
  af_vec[is.na(af_vec)] <- caf_af[is.na(af_vec)]
}

prep <- prepare_markers(
  gt          = gt,
  pos         = pos,
  ped         = ped,
  af          = af_vec,
  status      = status,
  chrom       = chrom,
  max_markers = max_markers
)

if (length(prep) == 0L) {
  cat("chrom\tmarker\tpos\taf\tname\n", file = file.path(out_dir, "markers.tsv"))
  cat("n_branches\tn_markers\n0\t0\n", file = file.path(out_dir, "counts.tsv"))
  quit(save = "no", status = 0L)
}

gt_clean <- S4Vectors::mcols(prep)$gt
pos_sub  <- GenomicRanges::start(prep)
af_sub   <- as.numeric(S4Vectors::mcols(prep)$af)
mnames   <- as.character(S4Vectors::mcols(prep)$name)
n_m      <- length(prep)

# Write master markers.tsv
markers_df <- data.frame(
  chrom  = rep(chrom, n_m),
  marker = paste0("M", seq_len(n_m)),
  pos    = pos_sub,
  af     = round(af_sub, 6L),
  name   = mnames,
  stringsAsFactors = FALSE
)
utils::write.table(
  markers_df,
  file      = file.path(out_dir, "markers.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

common_samples <- colnames(gt_clean)
seq_carriers   <- intersect(carriers, common_samples)
if (length(seq_carriers) == 0L) seq_carriers <- common_samples

branches <- split_pedigree(ped, targets = seq_carriers, max_bits = 22L)
valid_branches_count <- 0L

for (b_i in seq_along(branches)) {
  p_sub <- branches[[b_i]]
  if (!pedtools::is.ped(p_sub)) next
  b_ids <- as.character(labels(p_sub))
  if (length(intersect(b_ids, seq_carriers)) < 2L) next

  valid_m <- write_merlin_branch(
    p_sub      = p_sub,
    gt_clean   = gt_clean,
    pos_sub    = pos_sub,
    af_sub     = af_sub,
    carriers   = carriers,
    chrom      = chrom,
    out_dir    = out_dir,
    branch_idx = b_i,
    family_id  = famid_val
  )
  if (length(valid_m) >= 2L) {
    valid_branches_count <- valid_branches_count + 1L
  }
}

# Write summary counts
cat(sprintf("n_branches\tn_markers\n%d\t%d\n", valid_branches_count, n_m), file = file.path(out_dir, "counts.tsv"))
