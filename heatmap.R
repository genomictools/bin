#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort    <- args[1]
feature   <- args[2]
cnv       <- args[3]
type      <- args[4]
n_samples <- args[5]
n_genes   <- args[6]

# cnv   <- 'call-cnv-arrays/tests/results/scanned/pheno2.gene.cnv'
# density <- 0
# output  <- 'test.png'

# load cnv
d <- readr::read_table(cnv, col_names = FALSE)
d <- dplyr::select(d, region = X1, cn = X4, gene = X9, file = X5)
d <- dplyr::filter(d, gene != 'NOT_FOUND')
d <- dplyr::mutate(d, cn = as.integer(stringr::str_remove_all(cn, '.*=')))
d <- transform(d, gene = strsplit(gene, ','))
d <- tidyr::unnest(d, gene)
d <- tidyr::separate(d, file, c('cohort', 'key'), remove = FALSE)
d <- tidyr::separate(d, region, c('chrom', 'start', 'end'), remove = FALSE)

# columns: sample_count
sample_count <- table(d$file)
ca <- ComplexHeatmap::columnAnnotation(
  '# Samples' = ComplexHeatmap::anno_barplot(as.integer(sample_count))
)

# rows: gene_count
gene_count <- table(d$gene)
ra <- ComplexHeatmap::rowAnnotation(
  '# Variants' = ComplexHeatmap::anno_barplot(as.integer(gene_count))
)

# heatmap: cnvs
heat_matrix <- reshape2::acast(
  d,
  gene ~ key,
  value.var = 'cn',
  fill = '2',
  fun.aggregate = function(x) paste(unique(x), collapse = ',')
)

colors <- c("darkblue", "blue", "white", "red", "darkred") 
names(colors) <- c('0', '1', '2', '3', '4')

file_name <- paste(cohort, feature, type, 'png', sep = '.')

png(filename = file_name)

ComplexHeatmap::Heatmap(
  heat_matrix,
  col = colors,
  left_annotation = ra,
  top_annotation = ca,
  name = 'CNV',
  rect_gp = grid::gpar(col = "white", lwd = 2)
)

dev.off()
