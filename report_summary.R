#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort <- args[1]
type   <- args[2]
stats  <- args[3]

# read stats
stats <- readr::read_tsv(stats)

# plot by type
if ( type == 'samples' ) {
  # test samples
  d <- dplyr::filter(stats, selected)
  d <- dplyr::group_by(d, test.samples)
  d <- dplyr::summarise_if(d, is.numeric, mean)
  d <- dplyr::select(d, is.numeric)
  d <- as.matrix(d)
  
  g <- GGally::ggpairs(d)
  ggplot2::ggsave(
    plot = g,
    filename = paste(cohort, type, 'test', 'png', sep = '.')
  )
  
  # reference samples
  d <- dplyr::select(stats, -test.samples)
  d <- tidyr::pivot_longer(
      d,
      names_to = 'name',
      values_to = 'value',
      cols = c('correlations', 'expected.BF', 'phi', 'RatioSd', 'mean.p', 'median.depth')
    )
  g <- ggplot2::ggplot(d, ggplot2::
                         aes(y = ref.samples, x = value))
  g <- g + ggplot2::geom_point()
  g <- g + ggplot2::facet_grid(selected~name, scales = 'free_x')

  ggplot2::ggsave(
    plot = g,
    filename = paste(cohort, type, 'reference', 'png', sep = '.')
  )
} else if ( type == 'calls' ) {
  d <- dplyr::mutate(
    stats,
    width = end - start,
    type = ifelse(width > 1000, 'large', 'small')
  )
  for (l in c('small', 'large')) {
    d2 <- dplyr::filter(d, type == l)
    d2 <- dplyr::select(d2, nexons, width, BF, dplyr::starts_with('reads'))
    g <- GGally::ggpairs(d2)
    
    ggplot2::ggsave(
      plot = g,
      filename = paste(cohort, type, l, 'png', sep = '.')
    )
  }
}
