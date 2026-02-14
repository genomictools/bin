#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

gene   <- args[1]
key    <- args[2]
exons  <- args[3]
bam    <- args[4]
output <- args[5]

# load regions
exons <- readr::read_tsv(exons)
exons <- GenomicRanges::makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
exons <- GenomicRanges::reduce(exons)

# load bam and calculate coverage
bam <- unlist(strsplit(bam, ','))

test_coverage <- bamsignals::bamCoverage(bam[1], exons)
test_coverage <- unlist(as.list(test_coverage))

ref_coverage <- purrr::map(bam[-1], ~unlist(as.list(bamsignals::bamCoverage(.x, exons))))

# plot signal
ylim <- c(0, max(c(test_coverage, unlist(ref_coverage))))

# profile <- bamProfile(test, exons)
# counts <- bamCount(test, exons)

png(filename = output, width = 10, height = 5, units = 'in', res = 300)
plot(1:length(test_coverage),
     test_coverage,
     type = 'l',
     col = 'red',
     ylim = ylim,
     xlab = 'Relative Position (Exons)',
     ylab = 'Coverage',
     main = paste0(gene, ' (', key, ')'))

for (i in ref_coverage) {
  lines(1:length(test_coverage),
       i,
       type = 'l')
}

dev.off()
