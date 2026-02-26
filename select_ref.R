#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

test_key    <- args[1]
test_counts <- args[2]
ref_key     <- args[3]
ref_counts  <- args[4]
bins        <- args[5]
output      <- args[6]

# load tests
test <- rtracklayer::import.bed(test_counts)

# load ref
ref <- unlist(strsplit(ref_counts, ','))
ref <- purrr::map(ref, rtracklayer::import.bed)
ref <- purrr::map(ref, ~.x$score)
ref <- dplyr::bind_cols(ref)
ref <- as.matrix(ref)
colnames(ref) <- unlist(strsplit(ref_key, ','))
rownames(ref) <- test$name

# calculate correlation
optimal_ref <- ExomeDepth::select.reference.set(
  test.counts = test$score,
  reference.counts = ref,
  bin.length = GenomicRanges::width(test),
  n.bins.reduced = as.integer(bins)
)

# tidy stats
stats <- optimal_ref$summary.stats
stats <- dplyr::mutate(
  stats,
  selected = ifelse(ref.samples %in% optimal_ref$reference.choice, TRUE, FALSE),
  test.samples = unique(test_key)
)

stats <- dplyr::select(stats, test.samples, dplyr::everything())

readr::write_tsv(
  stats,
  output
)
