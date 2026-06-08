#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid    <- args[1]
rlist    <- args[2]
cases    <- args[3]
ped_file <- args[4]

# cases
cases <- readr::read_lines(cases)
# cases <- stringr::str_split(cases, '_', simplify = TRUE)[, 2]

# rlist
rlist <- unlist(strsplit(rlist, ','))
rlist <- purrr::map_df(rlist, ~readr::read_delim(.x, delim = ' ', col_names = c('variant', 'genotype', 'alt', 'ref')))
rlist <- tidyr::unite(rlist, samples, dplyr::starts_with('X'), sep = ' ')
rlist <- dplyr::mutate(rlist, samples = purrr::map_chr(stringr::str_split(samples, ' '), ~{paste(intersect(cases, unlist(.x)), collapse = ',')}))
rlist <- dplyr::select(rlist, variant, genotype, samples)
rlist <- transform(rlist, samples = strsplit(samples, ','))
rlist <- tidyr::unnest(rlist, samples)
rlist <- dplyr::mutate(rlist, genotype = ifelse(genotype == 'HET', 'a/b', 'a/a'))
rlist <- unique(rlist)

# Load pedigree
pdg  <- pedtools::readPed(ped_file, colSkip = c(6, 7))

rlist <- split(rlist, rlist$variant)

mms <- purrr::imap(rlist, ~{
    d <- tidyr::separate(
      .x, 
      variant,
      into = c('chrom', 'pos', 'ref', 'alt'),
      remove = FALSE
    )
    
    with(d, {
      s <- rep('a/a', length(cases))
      names(s) <- cases
      s[samples] <- genotype
      
      pedtools::marker(
        pdg,
        geno = s,
        name = unique(variant),
        chrom = unique(chrom),
        posMb = unique(as.integer(pos) / 1000000)
      )
    })
  })

pdg2 <- pedtools::setMarkers(pdg, mms)

out_file <- paste(famid, 'marked', sep = '.')
pedtools::writePed(pdg2, out_file)

aff  <- readr::read_delim(ped_file, col_select = c(2,6), col_names = FALSE)
aff <- dplyr::filter(aff, X6 == 2)
aff <- dplyr::pull(aff, X2)
readr::write_lines(aff, paste(famid, 'aff', 'txt', sep = '.'))

carr <- readr::read_delim(ped_file, col_select = c(2,7), col_names = FALSE)
carr <- dplyr::filter(carr, X7 == 3)
carr <- dplyr::pull(carr, X2)
readr::write_lines(carr, paste(famid, 'carr', 'txt', sep = '.'))

readr::write_lines(cases, paste(famid, 'star', 'txt', sep = '.'))
