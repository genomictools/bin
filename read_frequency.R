#!/usr/bin/env Rscript

# Capture command-line arguments
args        <- commandArgs(trailingOnly = TRUE)
type        <- args[1]
cohort      <- args[2]
category    <- args[3]
file        <- args[4]
size        <- args[5]
output_file <- args[6] 

# Functions
read_frequency <- function(fl, category, size, type = 'controls') {
  # Read file, add and group by category info
  d <- readr::read_tsv(fl)
  d <- dplyr::mutate(d, category = category)
  d <- dplyr::group_by(d, category, gene)
  
  # Get the case counts
  if ( type == 'controls' ) {
    # Estiamte case count from dominant and recessive models from summary data
    d <- dplyr::summarise(
      d,
      CONTROL_NVAR = sum(nvar),
      CONTROL_COUNT_DOM = sum(ac),
      CONTROL_COUNT_REC = round((round(ac/an)^2) * af) + nhom,
      CONTROL_SIZE = size
    )
  } else if ( type == 'cases' ) {
    # Case counts from genotype counts
      d <- dplyr::summarise(
        d,
        CASE_NVAR = sum(nvar),
        CASE_COUNT_DOM = het + hom,
        CASE_COUNT_REC = hom + ch,
        CASE_SIZE = size
      )
  }

  d <- dplyr::ungroup(d)
  d <- unique(d)
  return(d)
}

# Call read_frequency
counts <- read_frequency(
  file,
  category = category,
  size = size,
  type = type
)

# Write output
readr::write_tsv(counts, output_file)
