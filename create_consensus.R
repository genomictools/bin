#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cnv   <- args[1]
cn    <- as.integer(args[2])
ratio <- as.numeric(args[3])
stats <- as.integer(args[4])
output<- args[5]

# load library
library(GenomicRanges)

# load cnvs
cnv   <- unlist(strsplit(cnv, ','))
cnv <- purrr::map_df(cnv, readr::read_tsv)
cnv <- dplyr::rename(cnv, state = cn)

# filter cnvs
cnv <- dplyr::filter(cnv, state != as.integer(cn))
cnv <- dplyr::filter(cnv, ratio < as.numeric(ratio) | ratio > as.numeric(ratio) + 1 | is.na(ratio))
cnv <- dplyr::filter(cnv, stats > as.numeric(stats) | is.na(stats))

# convert to granges
gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
grl <- split(gr, gr$tool)

# find recurrent ranges
ro <- CNVRanger::populationRanges(grl)

# find overlap with original ranges
ol <- GenomicRanges::findOverlaps(ro, gr)

# summarize overlap data
d <- as.data.frame(GenomicRanges::mcols(gr)[S4Vectors::subjectHits(ol),])
d$range <- S4Vectors::queryHits(ol)

d <- dplyr::group_by(d, range)
d <- dplyr::left_join(
  dplyr::summarise_at(d, dplyr::vars(c('gene', 'key', 'tool')), ~paste(unique(.x), collapse = ',')),
  dplyr::summarise_at(d, dplyr::vars(c('gene', 'key', 'tool')), ~length(unique(.x))),
  by = 'range',
  suffix = c('', '_n')
)
d <- dplyr::select(d, -range)

GenomicRanges::mcols(ro) <- cbind(GenomicRanges::mcols(ro), d)
res <- as.data.frame(ro)

readr::write_tsv(res, output)

