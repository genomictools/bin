#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cnv    <- args[1]
ccn    <- args[2]
cratio <- args[3]
cstats <- args[4]
cwidth <- args[5]
style  <- args[6]
output <- args[7]

# load library
library(GenomicRanges)

# load cnvs
cnv   <- unlist(strsplit(cnv, ','))
cnv <- purrr::map(cnv, readr::read_tsv)
cnv <- purrr::map(cnv, ~GenomicRanges::makeGRangesFromDataFrame(.x, keep.extra.columns = TRUE))
cnv <- purrr::map(cnv, ~{x <- .x; GenomeInfoDb::seqlevelsStyle(x) <- style; x})
cnv <- purrr::map_df(cnv, tibble::as_tibble)
cnv <- dplyr::rename(cnv, state = cn)

# filter cnvs
cnv <- dplyr::filter(cnv, state != as.integer(ccn))
cnv <- dplyr::filter(cnv, ratio < (1 - as.numeric(cratio)) | ratio > (1 + as.numeric(cratio)) | is.na(ratio))
cnv <- dplyr::filter(cnv, stats > as.numeric(cstats) | is.na(stats))
cnv <- dplyr::filter(cnv, width < as.integer(cwidth))

# convert to granges
gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
grl <- split(gr, gr$key)

# find recurrent ranges
ro <- CNVRanger::populationRanges(grl, mode = 'RO')

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

