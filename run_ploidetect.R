# Rscript scripts/run_ploidetect_minimal.R \
#   --case POG965 \
#   --somatic-id P02866 \
#   --normal-id P02590 \
#   --somatic-bam /path/to/tumour.bam \
#   --normal-bam /path/to/normal.bam \
#   --genome-fa /path/to/hg38.fa \
#   --array-positions resources/snp_arrays/hg38/SNP_array_positions.txt \
#   --cyto resources/hg38/cytobands.txt \
#   --output-dir results \
#   --scripts-dir scripts \
#   --window-threshold 100000 \
#   --qual 10 \
#   --maxd 1000000


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# -----------------------------
# Helpers
# -----------------------------
run_cmd <- function(cmd, step_name) {
  message(sprintf("\n[%s] %s", step_name, cmd))
  status <- system(cmd)
  if (status != 0) stop(sprintf("Step '%s' failed with exit code %s", step_name, status))
}

run_cmd_to_dt <- function(cmd, step_name, col_types = NA) {
  message(sprintf("\n[%s] Parsing into memory", step_name))
  lines <- system(cmd, intern = TRUE)
  if (length(lines) == 0) stop(sprintf("Step '%s' produced no output", step_name))
  fread(text = paste(lines, collapse = "\n"), colClasses = col_types)
}

must_exist <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) label <- path
    stop(sprintf("Missing required file: %s", label))
  }
}

dir_for <- function(path) dirname(normalizePath(path, winslash = "/", mustWork = FALSE))

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# CLI
# -----------------------------
opt_list <- list(
  make_option("--case", type = "character", help = "Case ID"),
  make_option("--somatic-id", type = "character", help = "Somatic sample/library ID"),
  make_option("--normal-id", type = "character", help = "Normal sample/library ID"),
  make_option("--somatic-bam", type = "character", help = "Path to somatic BAM"),
  make_option("--normal-bam", type = "character", help = "Path to normal BAM"),
  make_option("--genome-fa", type = "character", help = "Reference genome FASTA"),
  make_option("--array-positions", type = "character", help = "2-col SNP positions (chr, pos)"),
  make_option("--cyto", type = "character", default = "", help = "Cytoband file path (optional)"),
  make_option("--output-dir", type = "character", help = "Output root dir"),
  make_option("--scripts-dir", type = "character", default = "scripts", help = "Directory containing helper scripts"),
  make_option("--window-threshold", type = "integer", default = 100000, help = "Window threshold for make_windows.py"),
  make_option("--qual", type = "integer", default = 10, help = "samtools depth -Q"),
  make_option("--maxd", type = "integer", default = 1000000, help = "samtools depth -m"),
  make_option("--threads", type = "integer", default = 8, help = "Reserved; passed where relevant"),
  make_option("--keep-tmp", action = "store_true", default = FALSE, help = "Keep intermediate files")
)

opt <- parse_args(OptionParser(option_list = opt_list))

required <- c("case","somatic_id","normal_id","somatic_bam","normal_bam","genome_fa","array_positions","output_dir")
missing_req <- required[!vapply(required, function(k) !is.null(opt[[k]]) && nzchar(opt[[k]]), logical(1))]
if (length(missing_req) > 0) {
  stop(sprintf("Missing required args: %s", paste(missing_req, collapse = ", ")))
}

must_exist(opt$somatic_bam, "--somatic-bam")
must_exist(opt$normal_bam, "--normal-bam")
must_exist(opt$genome_fa, "--genome-fa")
must_exist(paste0(opt$genome_fa, ".fai"), "FASTA index (.fai) for --genome-fa")
must_exist(opt$array_positions, "--array-positions")
if (nzchar(opt$cyto)) must_exist(opt$cyto, "--cyto")

# Required external scripts
make_windows_py <- file.path(opt$scripts_dir, "make_windows.py")
get_af_bash     <- file.path(opt$scripts_dir, "get_allele_freqs.bash")
merge_loh_r     <- file.path(opt$scripts_dir, "merge_loh.R")
prep_r          <- file.path(opt$scripts_dir, "prep_ploidetect2.R")
run_pd_r        <- file.path(opt$scripts_dir, "run_ploidetect2.R")
cnv_r           <- file.path(opt$scripts_dir, "ploidetect_copynumber.R")

for (p in c(make_windows_py, get_af_bash, merge_loh_r, prep_r, run_pd_r, cnv_r)) must_exist(p, p)

# -----------------------------
# Paths
# -----------------------------
pair_tag <- sprintf("%s_%s", opt$somatic_id, opt$normal_id)
pair_root <- file.path(opt$output_dir, opt$case, pair_tag)
scratch <- tempdir()  # Use temporary directory for scratch files only
scratch_loh <- file.path(scratch, sprintf("loh_%s", Sys.getpid()))

safe_mkdir(pair_root)
safe_mkdir(scratch_loh)

# Intermediate files preserved only when explicitly needed for external tools
windows_txt   <- file.path(scratch, sprintf("windows_%s.txt", Sys.getpid()))
loh_raw       <- file.path(scratch_loh, "loh_raw.txt")
merged_bed    <- file.path(scratch, sprintf("merged_%s.bed", Sys.getpid()))

segmented_rds <- file.path(pair_root, "segmented.RDS")
plots_pdf     <- file.path(pair_root, "plots.pdf")
models_txt    <- file.path(pair_root, "models.txt")
meta_rds      <- file.path(pair_root, "meta.RDS")
cna_txt       <- file.path(pair_root, "cna.txt")
cna_plots_pdf <- file.path(pair_root, "cna_plots.pdf")
cna_cond_txt  <- file.path(pair_root, "cna_condensed.txt")
metadat_rds   <- file.path(pair_root, "metadat.rds")

# -----------------------------
# 1) Build windows from normal depth (once per case/normal)
# -----------------------------
cmd_windows <- sprintf(
  "bash -lc %s",
  shQuote(sprintf(
    "samtools depth -Q %d -m %d %s | awk 'BEGIN{FS=OFS=\"\\t\"}{print $1,$2,$2+1,$3}' | python3 %s - %d | cut -f1,2,3 | bedtools sort -i stdin > %s",
    opt$qual, opt$maxd, shQuote(opt$normal_bam), shQuote(make_windows_py), opt$window_threshold, shQuote(windows_txt)
  ))
)
run_cmd(cmd_windows, "1_build_windows")
must_exist(windows_txt)
on.exit(try(unlink(windows_txt, force = TRUE), silent = TRUE), add = TRUE)

# ----------------------------------
# 2) Tumour binned depth from windows (in-memory)
# ----------------------------------
cmd_tumour <- sprintf(
  "bash -lc %s",
  shQuote(sprintf(
    "samtools depth -Q %d -m %d -a %s | awk 'BEGIN{FS=OFS=\"\\t\"}{print $1,$2,$2+1,$3}' | sort -k1,1 -k2,2n | bedtools map -b stdin -a %s -c 4 -o mean",
    opt$qual, opt$maxd, shQuote(opt$somatic_bam), shQuote(windows_txt)
  ))
)
tumour_dt <- run_cmd_to_dt(cmd_tumour, "2_tumour_binned_depth")
setnames(tumour_dt, c("chr", "start", "end", "depth_tumour"))

# ----------------------------------
# 3) Normal binned depth from windows (in-memory)
# ----------------------------------
cmd_normal <- sprintf(
  "bash -lc %s",
  shQuote(sprintf(
    "samtools depth -Q %d -m %d -a %s | awk 'BEGIN{FS=OFS=\"\\t\"}{print $1,$2,$2+1,$3}' | sort -k1,1 -k2,2n | bedtools map -b stdin -a %s -c 4 -o mean",
    opt$qual, opt$maxd, shQuote(opt$normal_bam), shQuote(windows_txt)
  ))
)
normal_dt <- run_cmd_to_dt(cmd_normal, "3_normal_binned_depth")
setnames(normal_dt, c("chr", "start", "end", "depth_normal"))
# -----------------------------
# 4) LOH/BAF from het SNPs
# -----------------------------
cmd_loh_raw <- sprintf(
  "bash -lc %s",
  shQuote(sprintf(
    "bash %s %s %s %s %s %s",
    shQuote(get_af_bash),
    shQuote(opt$normal_bam),
    shQuote(opt$somatic_bam),
    shQuote(opt$genome_fa),
    shQuote(opt$array_positions),
    shQuote(file.path(scratch, "loh_tmp"))
  ))
)
run_cmd(cmd_loh_raw, "4_loh_raw")
must_exist(loh_raw)
on.exit(try(unlink(scratch_loh, recursive = TRUE, force = TRUE), silent = TRUE), add = TRUE)

# ----------------------------------
# 4b) LOH/BAF binning (in-memory)
# ----------------------------------
message("\n[4b_loh_baf_binned] Computing BAF from SNP allele frequencies")
loh_raw_dt <- fread(loh_raw, col.names = c("chr", "pos", "ref", "ref_count", "alt", "alt_count"))
loh_raw_dt <- loh_raw_dt[(ref_count != 0 & alt_count != 0)]
loh_raw_dt[, baf := ref_count / (ref_count + alt_count)]
loh_dt <- loh_raw_dt[, .(chr, pos, baf)]
setkey(loh_dt, chr, pos)

# ----------------------------------
# 5) GC per window (in-memory)
# ----------------------------------
cmd_gc <- sprintf(
  "bash -lc %s",
  shQuote(sprintf(
    "bedtools nuc -fi %s -bed %s | cut -f1,2,3,5 | tail -n +2",
    shQuote(opt$genome_fa), shQuote(windows_txt)
  ))
)
gc_dt <- run_cmd_to_dt(cmd_gc, "5_gc_per_window")
setnames(gc_dt, c("chr", "start", "end", "gc"))

# ----------------------------------
# 6) Merge matrices in memory, then write merged_bed
# ----------------------------------
message("\n[6_merge_matrix] Merging depth, BAF, and GC into single matrix")
windows_dt <- fread(windows_txt, col.names = c("chr", "start", "end"))
setkey(windows_dt, chr, start, end)

# Merge tumour and normal depth
merged_dt <- tumour_dt[]
merged_dt[normal_dt, depth_normal := i.depth_normal, on = .(chr, start, end)]

# Map LOH SNPs to nearest window
setkey(loh_dt, chr, pos)
loh_binned <- windows_dt[loh_dt, baf := i.baf, on = .(chr, start <= pos, end > pos), mult = "first"]
merged_dt[loh_binned[!is.na(baf)], baf := i.baf, on = .(chr, start, end)]

# Add GC content
merged_dt[gc_dt, gc := i.gc, on = .(chr, start, end)]

# Remove 'chr' prefix for compatibility and write
merged_dt[, chr := gsub("^chr", "", chr)]
fwrite(merged_dt, merged_bed, sep = "\t", col.names = FALSE)
on.exit(try(unlink(merged_bed, force = TRUE), silent = TRUE), add = TRUE)
# -----------------------------
# 7) prep_ploidetect2.R
# -----------------------------
cmd_prep <- if (nzchar(opt$cyto)) {
  sprintf("Rscript %s -i %s -o %s -c %s",
          shQuote(prep_r), shQuote(merged_bed), shQuote(segmented_rds), shQuote(opt$cyto))
} else {
  sprintf("Rscript %s -i %s -o %s",
          shQuote(prep_r), shQuote(merged_bed), shQuote(segmented_rds))
}
run_cmd(cmd_prep, "7_presegment")
must_exist(segmented_rds)

# -----------------------------
# 8) run_ploidetect2.R
# -----------------------------
cmd_run_pd <- sprintf(
  "Rscript %s -i %s -o %s -p %s -r %s",
  shQuote(run_pd_r), shQuote(segmented_rds), shQuote(models_txt), shQuote(plots_pdf), shQuote(meta_rds)
)
run_cmd(cmd_run_pd, "8_run_ploidetect")
must_exist(models_txt)
must_exist(plots_pdf)
must_exist(meta_rds)

# -----------------------------
# 9) ploidetect_copynumber.R
# -----------------------------
cmd_cnv <- sprintf(
  "Rscript %s -i %s -m %s -p %s -o %s -d %s %s",
  shQuote(cnv_r),
  shQuote(segmented_rds),
  shQuote(models_txt),
  shQuote(cna_plots_pdf),
  shQuote(cna_txt),
  shQuote(metadat_rds),
  if (nzchar(opt$cyto)) paste("-c", shQuote(opt$cyto)) else ""
)
run_cmd(cmd_cnv, "9_copynumber")
must_exist(cna_txt)

message("\nDone.")
message(sprintf("Final outputs:\n- %s\n- %s\n- %s\n- %s", cna_txt, cna_plots_pdf, cna_cond_txt, metadat_rds))
message(sprintf("Important intermediate (segmented copy number): %s", segmented_rds))

message("\nDone.")
message(sprintf("Final outputs:\n- %s\n- %s\n- %s\n- %s", cna_txt, cna_plots_pdf, cna_cond_txt, metadat_rds))
message(sprintf("Important intermediate (segmented copy number): %s", segmented_rds))

# Cleanup temporary LOH directory
if (!opt$keep_tmp) {
  unlink(scratch_loh, recursive = TRUE, force = TRUE)
} else {
  message("Temporary files retained for debugging in: ", scratch_loh)
}