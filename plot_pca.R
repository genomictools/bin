#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

ref     <- args[1]
cohort  <- args[2]
mode    <- args[3]
scaled  <- args[4]
pop     <- args[5]
dims    <- args[6]

# Load population data
populations <- read.table(pop, header = FALSE)
names(populations) <- c('FID', 'IID', 'pop', 'superpop', 'type', 'assigned')

# Load scaled data
if ( mode == 'mds' ) {
    col_names <- c('FID', 'IID','SOL', paste0('D', 1:as.integer(dims)))
} else {
    col_names <- c('FID', 'IID', paste0('D', 1:as.integer(dims)))
}
scaled_dims <- readr::read_delim(scaled, col_names = col_names)
scaled_dims <- dplyr::left_join(scaled_dims, populations)
scaled_dims <- dplyr::mutate(scaled_dims, color = superpop)

# Plot
pca_plot <- (
    ggplot2::ggplot(scaled_dims) +
    ggplot2::geom_point(
        data = dplyr::filter(scaled_dims, type == 'references'),
        ggplot2::aes(x = D1, y = D2, color = as.factor(superpop)),
        size = .5, alpha = .5
    ) +
    ggplot2::geom_point(
        data = dplyr::filter(scaled_dims, type == 'cases'),
        ggplot2::aes(x = D1, y = D2),
        shape = 3, alpha = .5
    ) +
    ggplot2::labs(color = '') +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.ticks.length = ggplot2::unit(.15, "cm"),
        axis.ticks = ggplot2::element_line(color = 'lightgrey'),
        axis.title.y = ggplot2::element_text(vjust = 2)
    )
)

# Save plot
ggplot2::ggsave(
    plot = pca_plot,
    filename = paste(ref, cohort, mode, 'png', sep = '.'),
    width = 4.5, height = 3.5, units = 'in', dpi = 300
)
