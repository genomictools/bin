#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Meiotic Recombination Breakpoints & Crossover Tracks

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(GenomicRanges)
library(rtracklayer)

args         <- commandArgs(trailingOnly = TRUE)
n_args       <- length(args)
markers_file <- args[1L]
out_png      <- args[n_args]
out_bed      <- args[n_args - 1L]
chrom        <- as.character(args[n_args - 2L])
segment_file <- args[n_args - 3L]
mid_files    <- args[seq_len(n_args - 4L)[-1L]]
flow_files   <- grep("\\.flow$", mid_files, value = TRUE)

# Read markers
pos_all <- integer(0)
if (file.exists(markers_file) && file.info(markers_file)$size > 0L) {
  m_df <- utils::read.table(markers_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(m_df) > 0L && "pos" %in% colnames(m_df)) pos_all <- as.integer(m_df$pos)
}

# Read candidate interval if available
cand_int <- NULL
if (file.exists(segment_file) && file.info(segment_file)$size > 0L) {
  seg_df <- tryCatch(
    utils::read.table(segment_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(seg_df) && nrow(seg_df) > 0L && ncol(seg_df) >= 3L) {
    s_val <- as.numeric(seg_df[[2L]][1L])
    e_val <- as.numeric(seg_df[[3L]][1L])
    if (is.finite(s_val) && is.finite(e_val) && s_val < e_val) {
      cand_int <- c(as.integer(s_val), as.integer(e_val))
    }
  }
}

res_df <- parse_merlin_flow(
  flow_files         = flow_files,
  pos_all            = pos_all,
  chrom              = chrom,
  candidate_interval = cand_int
)

# Export BED table
if (nrow(res_df) > 0L) {
  con <- file(out_bed, open = "wt")
  writeLines("#chrom\tstart\tend\tbranch\tindividual\torigin\tfrom_allele\tto_allele\tcrossover_mb\tflank_type", con)
  utils::write.table(res_df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)
} else {
  file.create(out_bed)
}

# Render crossover tracks plot
plot_breakpoints(
  res_df             = res_df,
  pos_all            = pos_all,
  chrom              = chrom,
  candidate_interval = cand_int,
  out_png            = out_png
)
