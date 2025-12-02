#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

file    <- args[1]
pfb     <- args[2]
output  <- args[3]

# load calls
file <- readr::read_tsv(file)
file <- dplyr::select(file, `SNP Name` = Name, `Log R Ratio`, `B Allele Frequency` = `B Allele Freq`)

pfb <- readr::read_tsv(pfb)
pfb <- dplyr::select(pfb, `SNP Name` = SNP, Chromosome = Chr, Position)

res <- dplyr::inner_join(pfb, file)
res <- na.omit(res)

readr::write_tsv(res, output)
