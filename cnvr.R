#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort <- args[1]
tool   <- args[2]
cnv    <- args[3]
overlap<- args[4]
output <- args[5]

# load calls
cnv <- unlist(strsplit(cnv, ','))
names(cnv) <- unlist(strsplit(tool, ','))

col_names <- c('region', 'numsnp', 'length', 'state', 'sample', 'startsnp', 'endsnp', 'conf')
numeric_columns <- c("numsnp", "length", "state", "conf")

cnv <- purrr::imap(
    cnv,
    ~ {
        gr <- cnvr::read_cnv(
            .x,
            col_names = col_names,
            numeric_columns = numeric_columns
        )

        gr$tool <- .y
        gr
    }
)

# combine calls
cnv <- GenomicRanges::GRangesList(cnv)
cnv <- unlist(cnv)

# split by sample
cnv$sample <- stringr::str_remove_all(cnv$sample, '.data.*$')
grl <- GenomicRanges::GRangesList(split(cnv, cnv$sample))
grl <- GenomicRanges::sort(grl)

# get recurrence
cnvr <- CNVRanger::populationRanges(
    grl,
    mode = 'RO',
    ro.thresh = as.integer(overlap)
)

# assign sample and tool info to recurrent segments
ol <- GenomicRanges::findOverlaps(cnvr, cnv)
ol_subject <- cnv[S4Vectors::subjectHits(ol)]
ol_query <- S4Vectors::queryHits(ol)

cnvr$sample <- purrr::map_chr(split(ol_subject$sample, ol_query), ~paste(unique(.x), collapse = ','))
cnvr$freq_sample <- purrr::map_int(split(ol_subject$sample, ol_query), ~length(unique(.x)))
cnvr$tool <- purrr::map_chr(split(ol_subject$tool, ol_query), ~paste(unique(.x), collapse = ','))
cnvr$freq_tool <- purrr::map_int(split(ol_subject$tool, ol_query), ~length(unique(.x)))

# write to file
cnvr <- tibble::as_tibble(cnvr)
readr::write_tsv(cnvr, output)
