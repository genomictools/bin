#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

famid    	<- args[1]
ped_file	<- args[2]
affected    <- args[3]
carrier	    <- args[4]
starred	    <- args[5]
gene 	    <- args[6]
variant 	<- args[7]

# Load data
pedigree <- pedtools::readPed(ped_file)
affected <- readr::read_lines(affected)
carrier  <- readr::read_lines(carrier)
starred  <- readr::read_lines(starred)
size <- pedtools::pedsize(pedigree)
size <- max(4, size / 1.5)
file_name <- paste(famid, gene, "png", sep = '.')

# Extract marker names
markers <- unlist(strsplit(variant, split = ','))
title <- paste(c(famid, gene, markers), collapse = '\n')

# Plot the pedigree  
png(file_name,
    width = size, height = size,
    units = 'in', res = 300)

plot(
	pedigree,
    aff = affected,
    carrier = carrier,
    starred = starred,
    marker = markers,
    margins = c(0.6, 1, length(markers) + 4, 1),
    title = title
)

dev.off()
