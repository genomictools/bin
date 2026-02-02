#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort <- args[1]
gene   <- args[2]
key    <- args[3]
counts <- args[4]
ref    <- args[5]
prob   <- args[6]

# load tests
counts <- readr::read_tsv(counts)
counts <- GenomicRanges::makeGRangesFromDataFrame(counts, keep.extra.columns = TRUE)

# load ref
ref <- unlist(strsplit(ref, ','))
ref <- purrr::map_df(ref, readr::read_tsv)
ref <- GenomicRanges::makeGRangesFromDataFrame(ref, keep.extra.columns = TRUE)

# reshape
ref <- tidyr::pivot_wider(
  as.data.frame(GenomicRanges::mcols(ref)),
  names_from = 'key',
  values_from = 'counts'
)

ref_matrix <- as.matrix(dplyr::select(ref, -dplyr::starts_with('exon'), -dplyr::starts_with('GC')))

# aggregate references
ref_agg <- apply(X = ref_matrix, MAR = 1, FUN = sum)

# covariants
covar <- data.frame(GC = counts$GC)

# create ExomeDepth object
library(ExomeDepth)
all_exons <- new(
  'ExomeDepth',
  data = covar,
  test = counts$counts,
  reference = ref_agg,
  formula = 'cbind(test, reference) ~ GC'
)

# call CNVs using HMM
cnv <- ExomeDepth::CallCNVs(
  x = all_exons,
  transition.probability = as.numeric(prob),
  chromosome = as.character(GenomicRanges::seqnames(counts)),
  start = GenomicRanges::start(counts),
  end = GenomicRanges::end(counts),
  name = paste(gene, '-0_', 1:length(counts))
)

# save cnv object
readr::write_rds(
  cnv,
  paste(cohort, gene, key, 'cnv', 'rds', sep = '.')
)

# write output to text file
d <- tibble::as_tibble(cnv@CNV.calls)
if ( nrow(d) > 0) {
  d <- dplyr::mutate(d, cohort = cohort, gene = gene, key = key)
}

readr::write_tsv(
  d,
  paste(cohort, gene, key, 'cnv', 'tsv', sep = '.')
)
