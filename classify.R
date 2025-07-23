#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid       <- args[1]
category    <- args[2]
sharing     <- args[3]
type        <- args[4]

# read in the data
files <- unlist(strsplit(sharing, ','))
sharing <- purrr::map_df(files, readr::read_tsv, col_types = "cccccccccciiiiiiii")

# classify
res <- dplyr::mutate(sharing, sharing = dplyr::case_when(
    mac_non == 0 & expected_affected == 0 & (mac_obligate + mac_potential >= 1) ~ 'carriers',
    mac_non == 0 & mac_affected == 1 & (expected_affected + expected_obligate + expected_potential == 1) ~ 'singleton',
    mac_non == 0 & mac_affected == 1 & mac_obligate == 0 & mac_potential == 0 & mac_non == 0 ~ 'denovo',
    mac_non == 0 & mac_affected >= expected_affected & mac_obligate == expected_obligate ~ 'complete',
    mac_non == 0 & mac_affected < expected_affected & mac_obligate < expected_obligate ~ 'partial',
    TRUE ~ 'other'
  )) 

if ( type == 'gene' ) {
  res <- dplyr::filter(res, sharing != 'other')
  res <- dplyr::group_by(res, famid, category, ensembl, gene)
  res <- dplyr::reframe(
    res,
    nvar    = length(unique(variant)),
    nfam    = length(unique(famid)),
    variant = paste(unique(variant), collapse = ','),
    famid   = paste(unique(famid), collapse = ','),
    sharing = paste(unique(sharing), collapse = ',')
  )
}

res <- dplyr::mutate(res, type = type)

# Write output
output <- paste(famid, category, type, 'tsv', sep = '.')
readr::write_tsv(res, output)
