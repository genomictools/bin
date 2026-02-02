#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

test_counts <- args[1]
ref_counts  <- args[2]
bins        <- args[3]
output      <- args[4]

# load tests
test <- readr::read_tsv(test_counts)
test <- GenomicRanges::makeGRangesFromDataFrame(test, keep.extra.columns = TRUE)
wdth <- GenomicRanges::width(test)

# load ref
ref <- unlist(strsplit(ref_counts, ','))
ref <- purrr::map_df(ref, readr::read_tsv)
ref <- GenomicRanges::makeGRangesFromDataFrame(ref, keep.extra.columns = TRUE)

ref <- dplyr::select(as.data.frame(GenomicRanges::mcols(ref)), exon_id, key, counts)

# reshape
ref <- tidyr::pivot_wider(
  ref,
  names_from = 'key',
  values_from = 'counts'
)

ref <- as.matrix(dplyr::select(ref, -dplyr::starts_with('exon_id')))

# calculate correlation
optimal_ref <- ExomeDepth::select.reference.set(
  test.counts = test$counts,
  reference.counts = ref,
  bin.length = wdth,
  n.bins.reduced = as.integer(bins)
)

# tidy stats
stats <- optimal_ref$summary.stats
stats <- dplyr::mutate(
  stats,
  selected = ifelse(ref.samples %in% optimal_ref$reference.choice, TRUE, FALSE),
  test.samples = unique(test$key)
)

stats <- dplyr::select(stats, test.samples, dplyr::everything())

readr::write_tsv(
  stats,
  output
)
