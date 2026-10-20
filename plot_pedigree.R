#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Plot Pedigree Structure with Affection & Carrier Labels

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(pedtools)

args       <- commandArgs(trailingOnly = TRUE)
ped_file   <- args[1L]
pheno_file <- args[2L]
typed_file <- args[3L]
shrink     <- args[4L]
title      <- args[5L]
out_png    <- args[6L]

ped <- tryCatch(pedtools::readPed(ped_file), error = function(e) NULL)
if (is.null(ped)) {
  ped_df <- utils::read.table(ped_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ped_df) <- tolower(colnames(ped_df))
  if (ncol(ped_df) >= 5L && !("id" %in% colnames(ped_df))) {
    colnames(ped_df)[1:5] <- c("famid", "id", "fid", "mid", "sex")
  }
  ped <- pedtools::ped(id = as.character(ped_df$id), fid = as.character(ped_df$fid), mid = as.character(ped_df$mid), sex = as.integer(ped_df$sex))
}

all_ids  <- as.character(labels(ped))
affected <- character(0)
if (!is.null(pheno_file) && nzchar(pheno_file) && file.exists(pheno_file) && file.info(pheno_file)$size > 0L) {
  ph_df <- utils::read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ph_df) <- tolower(colnames(ph_df))
  col_aff <- if ("affected" %in% colnames(ph_df)) "affected" else if ("aff" %in% colnames(ph_df)) "aff" else colnames(ph_df)[2L]
  col_id  <- if ("id" %in% colnames(ph_df)) "id" else colnames(ph_df)[1L]
  affected <- intersect(as.character(ph_df[[col_id]][ph_df[[col_aff]] == 2]), all_ids)
} else if ("aff" %in% colnames(ped_df)) {
  affected <- intersect(as.character(ped_df$id[ped_df$aff == 2]), all_ids)
}

status   <- infer_lineage(ped, affected = affected)
carriers <- names(status)[status == 3L]

typed <- NULL
if (!is.null(typed_file) && file.exists(typed_file)) {
  lines <- trimws(readLines(typed_file, warn = FALSE))
  raw_typed <- lines[nzchar(lines)]
  if (length(raw_typed) > 0L) {
    matches <- list(
      intersect(raw_typed, all_ids),
      intersect(paste0("X", raw_typed), all_ids),
      intersect(sub("^X", "", raw_typed), all_ids)
    )
    counts <- vapply(matches, length, integer(1))
    if (max(counts) > 0L) typed <- matches[[which.max(counts)]]
  }
}
if (is.null(typed)) {
  typed <- all_ids[grepl("^[XS]", all_ids)]
}

if (isTRUE(shrink == "true")) {
  ped <- shrink_pedigree(ped, status = status, typed = typed)
}

plot_family_pedigree(
  ped      = ped,
  affected = affected,
  carriers = carriers,
  typed    = typed,
  title    = title,
  shrink   = isTRUE(shrink == "true"),
  out_png  = out_png
)
