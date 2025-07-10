#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid    <- args[1]
category <- args[2]
rlist    <- args[3]
cases    <- args[4]
ped_file <- args[5]

# cases
cases <- readr::read_lines(cases)

# rlist
rlist <- readr::read_delim(rlist, delim = ' ', col_names = c('variant', 'genotype', 'alt', 'ref'))
rlist <- tidyr::unite(rlist, samples, dplyr::starts_with('X'), sep = ' ')
rlist <- dplyr::mutate(rlist, samples = purrr::map_chr(stringr::str_split(samples, ' '), ~{paste(intersect(cases, unlist(.x)), collapse = ',')}))
rlist <- dplyr::select(rlist, variant, genotype, samples)
rlist <- tidyr::pivot_wider(rlist, names_from = 'genotype', values_from = 'samples')
rlist <- tidyr::gather(rlist, 'genotype', 'sample', -variant)
rlist <- transform(rlist, sample = strsplit(sample, ','))
rlist <- tidyr::unnest(rlist, sample)
rlist <- dplyr::mutate(rlist, genotype = ifelse(genotype == 'HET', 'a/b', 'b/b'))
rlist <- na.omit(rlist)

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
      s <- genotype
      names(s) <- sample
      
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

out_file <- paste(famid, category, 'marked', sep = '.')
pedtools::writePed(pdg2, out_file)

aff  <- readr::read_tsv(ped_file, col_select = 6)
readr::write_lines(dplyr::pull(aff), paste(famid, 'aff', 'txt', sep = '.'))

carr <- readr::read_tsv(ped_file, col_select = 7)
readr::write_lines(dplyr::pull(carr), paste(famid, 'carr', 'txt', sep = '.'))
