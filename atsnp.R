#!/usr/bin/env Rscript

args   <- commandArgs(trailingOnly = TRUE)
vcf    <- args[1]
motifs <- args[2]
output <- args[3]

# Load motif data
if ( motifs == 'ENCODE' ) {
  data(encode_library, package = 'atSNP')
  motifs <- encode_motif
} else if ( motifs == 'JASPAR' ) {
  data(jaspar_library, package = 'atSNP')
  motifs <- jaspar_motif
}

# Load variants
vcf <- VariantAnnotation::readVcfAsVRanges(vcf)

# Formate variants
d <- tibble::as_tibble(as.data.frame(vcf))
d <- dplyr::select(d, chr = seqnames, snp = start, a1 = alt, a2 = ref)
d <- unique(d)
d <- dplyr::mutate(d, snpid = paste(chr, snp, a1, a2, sep = ':'))
d <- dplyr::select(d, chr, snp, a1, a2, snpid)
d <- na.omit(d)
d <- dplyr::filter(d, stringr::str_length(d$a1) == 1 & stringr::str_length(d$a2) == 1)

# Load snps
# expects: snpid a1 a2   chr      snp
readr::write_delim(
  d, 
  'snp_tbl.txt',
  delim = '\t'
)

snp_info <- atSNP::LoadSNPData(
  'snp_tbl.txt',
  snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38",
  genome.lib = "BSgenome.Hsapiens.UCSC.hg38"
)

# Compute scores
atsnp.scores <- atSNP::ComputeMotifScore(motifs, snp_info, ncores = 1)

# Format scores
atsnp.scores <- dplyr::left_join(atsnp.scores$snp.tbl, atsnp.scores$motif.scores)
atsnp.scores <- tidyr::separate(atsnp.scores, snpid, into = c('CHROM','POS','REF','ALT'), sep = ':', remove = FALSE)
atsnp.scores <- dplyr::select(atsnp.scores, CHROM, POS, REF, ALT, ID = snpid, everything())

readr::write_tsv(atsnp.scores, output)
