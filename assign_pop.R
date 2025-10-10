#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

ref     <- args[1]
cohort  <- args[2]
mode    <- args[3]
scaled  <- args[4]
pop     <- args[5]
dims    <- args[6]

# ref     <- '1kg'
# cohort  <- 'fact'
# mode    <- 'noclusters'
# scaled  <- '../childhood_cancer/fact_ancestry/results/scaled/1kg.fact.noclusters.txt'
# pop     <- '../childhood_cancer/fact_ancestry/results/scaled/1kg.fact.pop'
# dims    <- 2

# Load population data
populations <- read.table(pop, header = FALSE)
names(populations) <- c('FID', 'IID', 'pop', 'superpop', 'type')
populations <- dplyr::mutate(populations, FID = as.character(FID))

# Load scaled data
if ( mode == 'mds' ) {
  col_names <- c('FID', 'IID','SOL', paste0('D', 1:as.integer(dims)))
} else {
  col_names <- c('FID', 'IID', paste0('D', 1:as.integer(dims)))
}

scaled_dims <- readr::read_delim(scaled, col_names = col_names)
scaled_dims <- dplyr::mutate(scaled_dims, FID = as.character(FID))
scaled_dims <- dplyr::left_join(scaled_dims, populations)

# Get known and unknown population centers
centers <- dplyr::mutate(scaled_dims, ID = ifelse(type == 'cases', IID, superpop))
centers <- dplyr::group_by(centers, ID)
centers <- dplyr::summarise_at(centers, dplyr::vars(dplyr::starts_with('D')), mean)
centers <- dplyr::ungroup(centers)

centers_mat <- as.matrix(dplyr::select(centers, dplyr::starts_with('D')))
rownames(centers_mat) <- centers$ID

distances <- as.matrix(dist(centers_mat))

ind <- rownames(distances) %in% na.omit(unique(scaled_dims$superpop))
distances <- distances[!ind, ind]

res <- tibble::tibble(
  IID = rownames(distances),
  assigned = colnames(distances)[apply(distances, 1, which.min)]
)

assigned_pop <- dplyr::left_join(populations, res)

readr::write_tsv(
  assigned_pop, 
  paste(ref, cohort, mode, 'pop', sep = '.'),
  col_names = FALSE
)
