#!/usr/bin/env Rscript
# Author: MOPOPGEN
# Task: Plot Theoretical vs Empirical Kinship Heatmap

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "_init.R"))
library(pedtools)

args        <- commandArgs(trailingOnly = TRUE)
genome_file <- args[1L]
ped_file    <- args[2L]
png_file    <- args[3L]

ped <- tryCatch(pedtools::readPed(ped_file), error = function(e) NULL)
if (is.null(ped)) {
  ped_df <- utils::read.table(ped_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ped_df) <- tolower(colnames(ped_df))
  if (ncol(ped_df) >= 5L && !("id" %in% colnames(ped_df))) {
    colnames(ped_df)[1:5] <- c("famid", "id", "fid", "mid", "sex")
  }
  ped <- pedtools::ped(id = as.character(ped_df$id), fid = as.character(ped_df$fid), mid = as.character(ped_df$mid), sex = as.integer(ped_df$sex))
}

genome <- tryCatch(
  utils::read.table(genome_file, header = TRUE, stringsAsFactors = FALSE),
  error = function(e) data.frame()
)

mkin_comb <- combine_kinship(genome, ped)
plot_kinship_heatmap(mkin_comb, out_png = png_file)
