#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
annotations <- args[1]
rlist       <- args[2]
output_file <- args[3]

# Load annotations
annotations <- readr::read_tsv(
  annotations,
  col_select = 1:2
)

# Load rlist
rlist <- readr::read_delim(
  rlist,
  delim = ' ',
  col_names = c('variant', 'genotype', 'ref', 'alt')
)

# Tidy up rlist
d <- tidyr::unite(rlist, 'samples', dplyr::starts_with('X'), sep = ' ')
d <- transform(d, samples = strsplit(samples, ' '))
d <- tidyr::unnest(d, samples)
d <- dplyr::filter(d, samples != '0')
d <- dplyr::filter(d, !is.na(samples))

# Join with annotations
d <- dplyr::inner_join(annotations, d)

# Number of variants
nvar <- dplyr::group_by(d, gene)
nvar <- dplyr::summarise(nvar, nvar = length(unique(variant)))

# Genotypes: het, hom, ch
het <- dplyr::filter(d, genotype == 'HET')
het <- dplyr::group_by(d, gene)
het <- dplyr::summarise(het, het = length(unique(samples)))

hom <- dplyr::filter(d, genotype == 'HOM')
hom <- dplyr::group_by(hom, gene)
hom <- dplyr::summarise(hom, hom = length(unique(samples)))

ch <- dplyr::group_by(d, gene, samples)
ch <- dplyr::summarise(ch, ch = length(unique(variant)))
ch <- dplyr::group_by(ch, gene)
ch <- dplyr::summarise(ch, ch = sum(ch > 1))

# Merge
res <- dplyr::left_join(nvar, het)
res <- dplyr::left_join(res, hom)
res <- dplyr::left_join(res, ch)
res <- dplyr::mutate_if(res, is.integer, ~ifelse(is.na(.x), 0, .x))

# Write to file
readr::write_tsv(
  res,
  output_file
)
