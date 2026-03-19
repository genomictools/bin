#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort <- args[1]
gene   <- args[2]
key    <- args[3]
tool   <- args[4]
test   <- args[5]
ref    <- args[6]
prob   <- args[7]

# load tests
test <- rtracklayer::import.bed(test)
names <- stringr::str_split(test$name, '\\-', simplify = TRUE)[, 3]
names <- paste0(
  names,
  '-0_',
  unname(unlist(purrr::map(table(names), ~seq(1, .x)))))

# load ref
ref <- unlist(strsplit(ref, ','))
ref <- purrr::map(ref, rtracklayer::import.bed)
ref <- purrr::map(ref, ~.x$score)
ref <- dplyr::bind_cols(ref)
ref <- as.matrix(ref)

# aggregate references
ref_agg <- apply(X = ref, MAR = 1, FUN = sum)

# create ExomeDepth object
library(ExomeDepth)
all_exons <- new(
  'ExomeDepth',
  test = test$score,
  reference = ref_agg,
  formula = 'cbind(test, reference) ~ 1'
)

# call CNVs using HMM
cnv <- ExomeDepth::CallCNVs(
  x = all_exons,
  transition.probability = as.numeric(prob),
  chromosome = as.character(GenomicRanges::seqnames(test)),
  start = GenomicRanges::start(test),
  end = GenomicRanges::end(test),
  name = names
)

# save cnv object
readr::write_rds(
  cnv,
  paste(cohort, gene, key, tool, 'cnv', 'rds', sep = '.')
)

# write output to text file
d <- tibble::as_tibble(matrix(ncol = 14))
colnames(d) <- c(
  'seqnames','start','end',
  'exon_start','exon_end','cn','cov','norm_cov','ratio','stats',
  'cohort','key','gene','tool'
)

if (nrow(cnv@CNV.calls) > 0) {
  d <- with(
    tibble::as_tibble(cnv@CNV.calls),
    tibble::tibble(
      seqnames = chromosome,
      start,
      end,
      exon_start = start.p,
      exon_end = end.p,
      cn = type,
      cov = reads.observed,
      norm_cov = reads.expected,
      ratio = reads.ratio,
      stats = BF,
      cohort = cohort,
      key = key,
      gene = gene,
      tool = tool
    )
  )
  d <- dplyr::mutate(d, cn = ifelse(cn == 'deletion', 1, 3))
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(cnv@CNV.calls)
  ol <- GenomicRanges::findOverlaps(gr, test)
  names <- stringr::str_split(test$name, '\\-', simplify = TRUE)[, 3]
  names <- ifelse(names == 'NA', NA, names)
  names <- split(
    names[S4Vectors::subjectHits(ol)],
    S4Vectors::queryHits(ol)
  )
  
  d$gene <- as.character(purrr::map(names, ~paste(na.omit(unique(.x)), collapse = ',')))
  
} else {
  d <- na.omit(d)
}

readr::write_tsv(
  d,
  paste(cohort, gene, key, tool, 'cnv', 'tsv', sep = '.')
)

# plot
for (g in unique(d$gene)) {
  anno <- dplyr::filter(cnv@annotations, stringr::str_detect(name, g))
  anno <- c(min(anno$start), max(anno$end))
  flank <- (anno[2] - anno[1])/2
  xlim <- c(anno[1] - flank, anno[2] + flank)

  png(paste(cohort, g, key, tool, 'cnv', 'png', sep = '.'))
  seg_cnv <- dplyr::filter(d, gene == g)
  ExomeDepth::plot(
    cnv,
    sequence = unique(as.character(seg_cnv$seqnames)),
    xlim = xlim,
    main = paste0(g, ' (', key, ')'),
    cex.lab = 0.8,
    with.gene = TRUE,
    xlab = unique(as.character(seg_cnv$seqnames)),
    count.threshold = 0
  )
  
  dev.off()
}
