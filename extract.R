#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

file     <- args[1]
name_col <- args[2]
baf_col  <- unlist(strsplit(args[3], split = ","))
lrr_col  <- unlist(strsplit(args[4], split = ","))
a1_col   <- unlist(strsplit(args[5], split = ","))
a2_col   <- unlist(strsplit(args[6], split = ","))
output   <- args[7]

# Get number of lines to skip
skip <- grep("\\[Data\\]", readr::read_lines(file, n_max = 100))
skip <- ifelse(length(skip) == 0, 0, skip)

# Read input file
res <- readr::read_tsv(file, skip = skip)

# Create GType column
if (a1_col != 'null' & a2_col != 'null') {
  res <- tidyr::unite(
    res,
    'GType',
    dplyr::any_of(c(a1_col, a2_col)),
    sep = ''
  )
}

# Subset columns
res <- dplyr::select(
  res,
  dplyr::any_of(c(name_col, baf_col, lrr_col, 'GType'))
)

# Rename columns
col_names <- c("Name", "B Allele Freq", "Log R Ratio")
if (ncol(res) > 3) { col_names <- c(col_names, 'GType') }
res <- setNames(res, col_names)

# Write to file
readr::write_tsv(res, output)
