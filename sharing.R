#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid       <- args[1]
category    <- args[2]
snplist     <- args[3]
rlist       <- args[4]
frq         <- args[5]
annotations <- args[6]
cases       <- args[7]
ped_file    <- args[8]
blacklist   <- args[9]

# blacklisted variants
blacklist <- readr::read_lines(blacklist)

# annotations
anno <- readr::read_tsv(annotations)
anno <- dplyr::select(anno, variant, gene = SYMBOL, ensembl = Gene, IMPACT, Consequence, clinsig = CLIN_SIG)

# cases
cases <- readr::read_lines(cases)

# rlist
rlist <- readr::read_delim(rlist, delim = ' ', col_names = c('variant', 'genotype', 'alt', 'ref'))
rlist <- tidyr::unite(rlist, samples, dplyr::starts_with('X'), sep = ' ')
rlist <- dplyr::mutate(rlist, samples = purrr::map_chr(stringr::str_split(samples, ' '), ~{paste(intersect(cases, unlist(.x)), collapse = ',')}))
rlist <- dplyr::select(rlist, variant, genotype, samples)
rlist <- tidyr::pivot_wider(rlist, names_from = 'genotype', values_from = 'samples')

cols <- c('HET', 'HOM')

m <- as.data.frame(matrix(NA, ncol = length(cols)))
names(m) <- cols
rlist <- dplyr::left_join(rlist, m)

# pedigree
clusters <- tibble::tibble(
  V1 = as.integer(0:3),
  V2 = c('non', 'potential', 'affected', 'obligate')
)

# expected numbers in each cluster
carr <- readr::read_tsv(ped_file, col_select = c(2, 7))
carr <- dplyr::filter(carr, id %in% cases)
carr <- dplyr::select(carr, carr)
carr <- dplyr::left_join(carr, clusters, by = c('carr' = 'V1'))
carr <- dplyr::group_by(carr, cluster = V2)
carr <- dplyr::reframe(carr, expected = dplyr::n())

# mac in each cluster
frq <- read.table(frq, skip = 1)
frq <- setNames(frq, c('chrom', 'variant', 'cluster', 'alt', 'ref', 'maf', 'mac', 'nchromobs'))
frq <- dplyr::left_join(frq, clusters, by = c('cluster' = 'V1'))
frq <- dplyr::select(frq, variant, mac, cluster = V2)
frq <- tibble::as_tibble(frq)

# add info
info <- tibble::tibble(famid = famid, category = category, variant = unique(frq$variant))

# merge
res <- dplyr::full_join(info, frq)
res <- dplyr::left_join(res, carr)
res <- dplyr::left_join(res, rlist)
res <- dplyr::inner_join(res, anno)
res <- dplyr::filter(res, !variant %in% blacklist)
res <- tidyr::pivot_wider(res, values_from = c('mac', 'expected'), names_from = 'cluster')

# order columns
cols <- c(paste('expected', clusters$V2, sep = '_'),
          paste('mac', clusters$V2, sep = '_'))

m <- as.data.frame(matrix(NA, ncol = length(cols)))
names(m) <- cols

res <- dplyr::left_join(res, m)
res <- dplyr::relocate(res, setdiff(names(res), cols), sort(cols))

# fill NA
res <- dplyr::mutate_at(res, dplyr::vars(dplyr::starts_with('expected_')), ~ifelse(is.na(.x), 0, .x))
res <- dplyr::mutate_at(res, dplyr::vars(dplyr::starts_with('mac_')), ~ifelse(is.na(.x), 0, .x))

# Write output
output <- paste(famid, category, 'tsv', sep = '.')
readr::write_tsv(res, output)
