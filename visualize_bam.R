#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort   <- args[1]
gene     <- args[2]
test_key <- args[3]
test_cov <- args[4]
ref_key  <- args[5]
ref_cov  <- args[6]
output   <- args[7]

# load regions
test_cov <- rtracklayer::import.bed(test_cov)
l <- length(test_cov)

ref_cov <- unlist(strsplit(ref_cov, ','))
ref_cov <- purrr::map(ref_cov, rtracklayer::import.bed)
names(ref_cov) <- unlist(strsplit(ref_key, ','))

scores <- purrr::map(ref_cov, ~.x$score)

# plot signal
ylim <- c(0, max(unlist(scores), test_cov$score))

png(filename = output, width = 10, height = 5, units = 'in', res = 300)
plot(1:l,
     test_cov$score,
     type = 'l',
     col = 'red',
     ylim = ylim,
     xlab = 'Relative Position (Exons)',
     ylab = 'Coverage',
     main = paste0(gene, ' (', test_key, ')'))

for (i in scores) {
  lines(1:l,
       i,
       type = 'l')
}

dev.off()
