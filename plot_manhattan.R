#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort     <- args[1]
category   <- args[2]
test       <- args[3]
phenotype  <- args[4]
test_file  <- args[5]

# cohort     <- 'pheno'
# test       <- 'assoc'
# test_file  <- 'identify-associated-loci/tests/results/tests/pheno.Unfiltered.P1.assoc'

# Load data
d <- readr::read_table(test_file)
d <- dplyr::select(d, !dplyr::starts_with('X'))

x <- 1:nrow(d)
y <- -log10(d$P)
c <- ifelse(d$CHR%%2 == 1, 'black', 'darkgray')

# Create and save plot
png(paste(cohort, category, test, phenotype, 'png', sep = '.'),
    height = 5, width = 5, units = 'in', res = 300)

plot(x, y,
     pch = 19,
     col = c,
     xlab = "Chromosome",
     ylab = "P-value (-log_10)",
     xaxt = 'n')

axis(1, at = which(!duplicated(d$CHR)), labels = unique(d$CHR))

dev.off()