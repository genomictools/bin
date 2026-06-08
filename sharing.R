#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid       <- args[1]
category    <- args[2]
rlist       <- args[3]
annotations <- args[4]
cases       <- args[5]
ped_file    <- args[6]
blacklist   <- args[7]

# blacklisted variants
blacklist <- readr::read_lines(blacklist)

# annotations
anno <- readr::read_tsv(annotations)
anno <- dplyr::select(anno, variant, gene = SYMBOL, ensembl = Gene, IMPACT, Consequence, clinsig = CLIN_SIG)

# cases
cases <- readr::read_lines(cases)
# cases <- stringr::str_split(cases, '_', simplify = TRUE)[, 2]

# pedigree
clusters <- tibble::tibble(
  V1 = as.integer(0:3),
  V2 = c('non', 'potential', 'affected', 'obligate')
)

# expected numbers in each cluster
id_carr <- readr::read_delim(ped_file, col_select = c(2, 7), col_types = 'cc', col_names = FALSE)
id_carr <- setNames(id_carr, c('id', 'carr'))

carr <- dplyr::filter(id_carr, id %in% cases)
carr <- dplyr::select(carr, carr)
carr <- dplyr::left_join(carr, clusters, by = c('carr' = 'V1'))
carr <- dplyr::group_by(carr, cluster = V2)
carr <- dplyr::reframe(carr, expected = dplyr::n())
carr <- tidyr::pivot_wider(carr, names_from = 'cluster', values_from = 'expected', names_prefix = 'expected_')
carr <- dplyr::mutate(carr, famid = famid)

# rlist
rlist <- readr::read_delim(rlist, delim = ' ', col_names = c('variant', 'genotype', 'alt', 'ref'))
rlist <- tidyr::unite(rlist, samples, dplyr::starts_with('X'), sep = ' ')
rlist <- dplyr::mutate(rlist, samples = purrr::map_chr(stringr::str_split(samples, ' '), ~{paste(intersect(cases, unlist(.x)), collapse = ',')}))
rlist <- dplyr::select(rlist, variant, genotype, samples)
rlist <- transform(rlist, samples = strsplit(samples, ','))
rlist <- tidyr::unnest(rlist, samples)

genotypes <- tidyr::pivot_wider(rlist, names_from = 'genotype', values_from = 'samples', values_fn = function(x) paste(x, collapse = ','))
cols <- c('HET', 'HOM')

m <- as.data.frame(matrix(NA, ncol = length(cols)))
names(m) <- cols
genotypes <- dplyr::left_join(genotypes, m)
genotypes <- dplyr::mutate_all(genotypes, ~ifelse(is.na(.x), '', .x))

# mac in each cluster
frq <- dplyr::left_join(rlist, id_carr, by = c('samples'='id'))
frq <- dplyr::left_join(frq, clusters, by = c('carr'='V1'))
frq <- dplyr::group_by(frq, variant, cluster = V2)
frq <- dplyr::reframe(frq, mac = length(unique(samples)))
frq <- tidyr::pivot_wider(frq, names_from = 'cluster', values_from = 'mac', names_prefix = 'mac_', values_fill = 0)

# add info
info <- tibble::tibble(famid = famid, category = category, variant = unique(frq$variant))

# merge
res <- dplyr::inner_join(info, anno)
res <- dplyr::filter(res, !variant %in% blacklist)
res <- dplyr::left_join(res, genotypes)
res <- dplyr::left_join(res, carr)
res <- dplyr::left_join(res, frq)

# order columns
cols <- c(paste('expected', clusters$V2, sep = '_'),
          paste('mac', clusters$V2, sep = '_'))

m <- as.data.frame(matrix(NA, ncol = length(cols)))
names(m) <- cols

res <- dplyr::left_join(res, m)
res <- dplyr::relocate(res, setdiff(names(res), cols), sort(cols))
res <- dplyr::mutate_all(res, ~ifelse(is.na(.x), 0, .x))

# Write output
output <- paste(famid, category, 'tsv', sep = '.')
readr::write_tsv(res, output)
