# Auto-load latest package development code if available
suppressPackageStartupMessages(library(linker))

# Resolve directory of the current script or _init.R
.get_script_dir <- function() {
  # Check --file argument
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa) > 0L) {
    p <- sub("^--file=", "", fa[1L])
    if (nzchar(p)) {
      ap <- normalizePath(p, mustWork = FALSE)
      if (file.exists(ap)) return(dirname(ap))
      # If relative to PATH
      sp <- Sys.which(basename(p))
      if (nzchar(sp)) return(dirname(normalizePath(sp, mustWork = FALSE)))
    }
  }
  # Fallback to sys.frame or default repository structure
  def_p <- "/data/rds/DGE/DUDGE/MOPOPGEN/mahmed03/genomictools/identify-family-haplotypes/bin"
  if (dir.exists(def_p)) return(def_p)
  getwd()
}

.sdir <- .get_script_dir()
.rdir <- normalizePath(file.path(.sdir, "..", "lib", "linker", "R"), mustWork = FALSE)
if (dir.exists(.rdir)) {
  for (.f in list.files(.rdir, "\\.R$", full.names = TRUE)) {
    source(.f)
  }
}

