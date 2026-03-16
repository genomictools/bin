#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cohort  <- args[1]
gene    <- args[2]
test_key<- args[3]
test    <- args[4]
ref_key <- args[5]
ref     <- args[6]
gc      <- args[7]
pmap    <- args[8]

# load tests
test <- rtracklayer::import.bed(test)
gc <- readr::read_csv(gc)
pmap <- readr::read_csv(pmap)
test$gc <- gc$score
test$mapp <- pmap$score

# load ref
ref <- unlist(strsplit(ref, ','))
ref <- purrr::map(ref, rtracklayer::import.bed)
ref <- purrr::map(ref, ~.x$score)
ref <- dplyr::bind_cols(ref)
ref <- as.matrix(ref)
ref_keys <- unlist(strsplit(ref_key, ','))

counts <- cbind(test$score, ref) + 1
colnames(counts) <- c(test_key, ref_keys)
rownames(counts) <- test$name

test <- test[rowSums(counts) > 10]
counts <- counts[rowSums(counts) > 10,]

# normalization
norm_index <- which(colnames(counts) %in% ref_keys)
chr_indices <- split(seq_len(nrow(counts)), as.character(GenomeInfoDb::seqnames(test)))
cnv_list <- list()

for (chr_name in names(chr_indices)) {
  idx <- chr_indices[[chr_name]]
  counts_chr <- counts[idx, , drop = FALSE]
  test_chr <- test[idx]

  if (nrow(counts_chr) < 2) {
    next
  }

  nonzero <- counts_chr[apply(counts_chr, 1, function(x) !any(x == 0)), , drop = FALSE]
  if (nrow(nonzero) == 0) {
    next
  }

  pseudo <- exp(rowMeans(log(nonzero)))
  N <- apply(sweep(nonzero, 1, pseudo, "/"), 2, median)

  K_grid <- 1:min(5, ncol(counts_chr) - 1)
  if (length(K_grid) == 0) {
    next
  }
  
  set.seed(1234)
  norm_ns <- CODEX2::normalize_codex2_ns(
    Y_qc = counts_chr,
    gc_qc = test_chr$gc,
    K = K_grid,
    norm_index = norm_index,
    N = N
  )

  optK <- K_grid[which.max(norm_ns$BIC)]

  set.seed(1234)
  seg_cbs <- CODEX2::segmentCBS(
    counts_chr[, test_key, drop = FALSE],
    norm_ns$Yhat,
    optK = optK,
    K = K_grid,
    sampname_qc = colnames(counts_chr),
    ref_qc = IRanges::ranges(test_chr),
    chr = as.character(GenomeInfoDb::seqnames(test_chr)),
    lmax = 400,
    mode = "integer"
  )

  cnv_chr <- tibble::as_tibble(seg_cbs)
  if (nrow(cnv_chr) == 0) {
    next
  }

  cnv_list[[chr_name]] <- with(
    cnv_chr,
    tibble::tibble(
      seqnames = chr,
      start = st_bp,
      end = ed_bp,
      exon_start = st_exon,
      exon_end = ed_exon,
      cn = copy_no,
      cov = raw_cov,
      norm_cov = norm_cov,
      ratio = lratio,
      stats = mBIC,
      cohort = cohort,
      key = test_key,
      gene = gene
    )
  )
}

if (length(cnv_list) > 0) {
  cnv <- dplyr::bind_rows(cnv_list)
  cnv <- stats::na.omit(cnv)
} else {
  cnv <- tibble::tibble(
    seqnames = character(),
    start = integer(),
    end = integer(),
    exon_start = integer(),
    exon_end = integer(),
    cn = numeric(),
    cov = numeric(),
    norm_cov = numeric(),
    ratio = numeric(),
    stats = numeric(),
    cohort = character(),
    key = character(),
    gene = character()
  )
}
cnv <- dplyr::filter(cnv, end - start > 1)
gr <- GenomicRanges::makeGRangesFromDataFrame(cnv)
ol <- GenomicRanges::findOverlaps(gr, test)
names <- stringr::str_split(test$name, '\\-', simplify = TRUE)[, 3]
names <- ifelse(names == 'NA', NA, names)
names <- split(
  names[S4Vectors::subjectHits(ol)],
  S4Vectors::queryHits(ol)
)
cnv$gene <- as.character(purrr::map(names, ~na.omit(unique(.x))))

readr::write_tsv(
  cnv,
  paste(cohort, gene, test_key, 'cnv', 'tsv', sep = '.')
)

# plot
for( g in unique(cnv$gene)) {
  seg_cnv <- dplyr::filter(cnv, gene == g)
  seg_counts <- counts[stringr::str_detect(test$name, g),]

  png(paste(cohort, g, test_key, 'cnv', 'png', sep = '.'))
  
  plot(
    1:nrow(seg_counts),
    ylim = c(min(seg_counts), max(seg_counts)),
    type = 'n',
    xlab = 'Exome number',
    ylab = 'Read Count',
    main = paste0(g, ' (', test_key, ')')
  )
  
  for ( i in unlist(strsplit(ref_key, ',')) ) {
    lines(
      1:nrow(seg_counts),
      seg_counts[, i],
      col = 'grey'
    )
  }
  
  lines(
    1:nrow(seg_counts),
    seg_counts[, test_key]
  )
  
  points(
    1:nrow(seg_counts),
    seg_counts[, test_key],
    pch = 19,
    cex = 0.5
  )
  
  for (i in seq_len(nrow(seg_cnv))) {
    with(seg_cnv[i,],
         {
           rect(
             exon_start, 0, exon_end, max(seg_counts),
             border = FALSE,
             col = adjustcolor(ifelse(cn == 3, 'red', 'green'), alpha.f = .2)
           )
         })
  }

  dev.off()
}
