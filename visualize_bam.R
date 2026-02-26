#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort   <- args[1]
gene     <- args[2]
key      <- args[3]
samples  <- args[4]
coverage <- args[5]
output   <- args[6]

# load regions
coverage <- unlist(strsplit(coverage, ','))
coverage <- purrr::map(coverage, rtracklayer::import.bed)
names(coverage) <- unlist(strsplit(samples, ','))

scores <- purrr::map(coverage, ~.x$score)
ind <- which(names(scores) == key)
l <- length(scores[[ind]])

# plot signal
ylim <- c(0, max(unlist(scores)))

png(filename = output, width = 10, height = 5, units = 'in', res = 300)
plot(1:l,
     scores[[ind]],
     type = 'l',
     col = 'red',
     ylim = ylim,
     xlab = 'Relative Position (Exons)',
     ylab = 'Coverage',
     main = paste0(gene, ' (', key, ')'))

for (i in scores[-ind]) {
  lines(1:l,
       i,
       type = 'l')
}

dev.off()
