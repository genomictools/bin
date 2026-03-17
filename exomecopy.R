#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(S4Vectors)
})

if (!exists('Rle', mode = 'function')) {
  stop('S4Vectors::Rle is not available in this runtime environment.')
}

cohort      <- args[1]
gene        <- args[2]
test_key    <- args[3]
tool        <- args[4]
test_counts <- args[5]
ref_key     <- args[6]
ref_counts  <- args[7]
gc          <- args[8]
pmap        <- args[9]
prob        <- args[10]

# load references
ref_counts <- purrr::map(unlist(strsplit(ref_counts, ',')), rtracklayer::import.bed)
ref_counts <- purrr::map(ref_counts, ~.x$score)
ref_counts <- dplyr::bind_cols(ref_counts)
names(ref_counts) <- unlist(strsplit(ref_key, ','))

# load counts
counts <- rtracklayer::import.bed(test_counts)
names(GenomicRanges::mcols(counts)) <- c('name', test_key)
GenomicRanges::mcols(counts) <- cbind(
  GenomicRanges::mcols(counts),
  ref_counts
)

gc <- rtracklayer::import.bed(gc)

# generate background
bg = as.numeric(exomeCopy::generateBackground(unlist(strsplit(ref_key, ',')), counts, median))

GenomicRanges::mcols(counts) <- cbind(
  GenomicRanges::mcols(counts),
  bg = bg,
  log.bg = log(bg + .1),
  width = GenomicRanges::width(counts),
  gc = gc$score
)

# call CNVs
library(S4Vectors)

cnv <- lapply(test_key, function(sample.name) {
  lapply(GenomeInfoDb::seqlevels(counts), function(seq.name) {
    exomeCopy::exomeCopy(
      GenomicRanges::sort(counts[GenomeInfoDb::seqnames(counts) == seq.name]),
      sample.name,
      X.names = c("log.bg", "width", "gc"),
      goto.cnv = as.numeric(prob)
    )
  })
})

# save cnv object
readr::write_rds(
  cnv,
  paste(cohort, gene, test_key, tool,  'cnv', 'rds', sep = '.')
)

# write output to text file
cnv <- exomeCopy::compileCopyCountSegments(cnv)

ol <- GenomicRanges::findOverlaps(cnv, counts)
names <- stringr::str_split(counts$name, '\\-', simplify = TRUE)[, 3]
names <- ifelse(names == 'NA', NA, names)

names <- split(
  names[S4Vectors::subjectHits(ol)],
  S4Vectors::queryHits(ol)
)
names <- as.character(purrr::map(names, ~paste(na.omit(unique(.x)), collapse = ',')))
cnv$gene <- names

d <- with(
  tibble::as_tibble(cnv),
  tibble::tibble(
    seqnames,
    start,
    end,
    exon_start = NA,
    exon_end = NA,
    cn = copy.count,
    cov = NA,
    norm_cov = NA,
    ratio = NA,
    stats = log.odds,
    cohort = cohort,
    key = sample.name,
    gene,
    tool = tool
  )
)

readr::write_tsv(
  d,
  paste(cohort, gene, test_key, tool, 'cnv', 'tsv', sep = '.')
)

# plot
purrr::map(
  GenomeInfoDb::seqlevels(cnv),
  ~{
    png(paste(cohort, .x, test_key, tool, 'cnv', 'png', sep = '.'))
    
    exomeCopy::plotCompiledCNV(
      cnv,
      seq.name = .x
    )
    
    dev.off()
  }
)
