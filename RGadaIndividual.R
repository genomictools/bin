#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
key             <- args[1]
input           <- args[2]
a_alpha         <- args[3]
t_statistic     <- args[4]
min_seg_length  <- args[5]
output          <- args[6]

# Get chrom names, and exclude if has a few markers
chrs <- read.delim(input)
chrs <- chrs$Chromosome
chrs <- names(table(chrs)[table(chrs) > as.integer(min_seg_length)])

# Setup
cnv.call <- gada::setupGADA(
  input,
  log2ratioCol = 4,
  BAFcol = 5,
  chrs = chrs
) 

# SBL
cnv.call <- gada::SBL(
  cnv.call,
  estim.sigma2 = TRUE,
  aAlpha = as.numeric(a_alpha),
  verbose = TRUE
)

# BE
cnv.call <- gada::BackwardElimination(
  cnv.call,
  T = t_statistic,
  MinSegLen = min_seg_length
)

# Summary
res <- summary( cnv.call )
# cnvs <- summary(
#   cnv.call,
#   length.base = c(500,10e6)
# )

# Tidy
res <- dplyr::mutate(
  tibble::as_tibble(res),
  sample = key,
  sample_index = sample,
  copy_number = State + 2,
  size = EndProbe - IniProbe,
  per_probe_score = abs(MeanAmp),
  lod_score = per_probe_score
)

res <- dplyr::select(
  res,
  sample,
  sample_index,
  copy_number,
  chr = chromosome,
  start = IniProbe,
  end = EndProbe,
  per_probe_score,
  size,
  num_probes = LenProbe,
  lod_score
)

# Write to files
write.table(
  res,
  file = output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
