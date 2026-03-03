#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort  <- args[1]
gene    <- args[2]
test_key<- args[3]
test    <- args[4]
ref_key <- args[5]
ref     <- args[6]

# load tests
test <- rtracklayer::import.bed(test)
test$gc <- rnorm(length(test), 0.5, .1)
test$mapp <- rnorm(length(test), 0.8, .1)

# load ref
ref <- unlist(strsplit(ref, ','))
ref <- purrr::map(ref, rtracklayer::import.bed)
ref <- purrr::map(ref, ~.x$score)
ref <- dplyr::bind_cols(ref)
ref <- as.matrix(ref)

counts <- cbind(test$score, ref)
colnames(counts) <- c(test_key, unlist(strsplit(ref_key, ',')))
rownames(counts) <- test$name

test <- test[rowSums(counts) > 10]
counts <- counts[rowSums(counts) > 10,]

# normalization
nonzero <- counts[apply(counts, 1, function(x) !any(x == 0)), ]
pseudo <- apply(nonzero, 1, function(x) exp(mean(log(x))))
N <- apply(apply(nonzero, 2, function(x) x / pseudo), 2, median)

norm_ns <- CODEX2::normalize_codex2_ns(
  Y_qc = counts,
  gc_qc = test$gc,
  K = 1:5,
  norm_index = which(colnames(counts) %in% unlist(strsplit(ref_key, ','))),
  N = N
)

# segmentation
seg_cbs <- CODEX2::segmentCBS(
  counts[, test_key, drop = FALSE],
  norm_ns$Yhat,
  optK = which.max(norm_ns$BIC),
  K = 1:5,
  sampname_qc = colnames(counts),
  ref_qc = IRanges::ranges(test),
  chr = GenomeInfoDb::seqlevels(test),
  lmax = 400,
  mode = "integer"
)

# write output to text file
readr::write_tsv(
  as.data.frame(seg_cbs),
  paste(cohort, gene, test_key, 'cnv', 'tsv', sep = '.')
)
