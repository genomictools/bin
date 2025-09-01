#!/usr/bin/env Rscript

library(gada)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), action = "store", type = "character", default = NA,
              help = "Input .baflrr file"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = "results",
              help = "Name of the output file with cnv calls"),
  make_option(c("-t", "--t_statistic"), action = "store", type = "numeric", default = 4,
              help = "Threshold for the t-statistic"),
  make_option(c("-a", "--a_alpha"), action = "store", type = "numeric", default = 0.8,
              help = "A alpha"),
  make_option(c("-m", "--min_seg_length"), action = "store", type = "numeric", default = 100,
              help = "Minimum segment length")
)

opt = parse_args(OptionParser(option_list = option_list))

cnv.call <- setupGADA(
  opt$input,
  log2ratioCol = 4,
  BAFcol = 5
) 

cnv.call <- SBL(
  cnv.call,
  estim.sigma2 = TRUE,
  aAlpha = opt$a_alpha,
  verbose = TRUE
)

cnv.call <- BackwardElimination(
  cnv.call,
  T = opt$t_statistic,
  MinSegLen = opt$min_seg_length
)

cnvs <- summary( cnv.call )
# cnvs <- summary(
#   cnv.call,
#   length.base = c(500,10e6)
# )

cnvs <- as_tibble(cnvs)
cnvs <- filter(cnvs, State != 0)
cnvs <- mutate(
  cnvs,
  sample = opt$input, sample_index = sample,
  copy_number = ifelse(State == 1, 3, 1),
  size = EndProbe - IniProbe,
  per_probe_score = abs(MeanAmp),
  lod_score = per_probe_score
)

d <- dplyr::select(
  cnvs,
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

write.table(
  d,
  file = opt$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
