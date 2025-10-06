#!/usr/bin/env Rscript

# Capture command-line arguments
args        <- commandArgs(trailingOnly = TRUE)
cases       <- args[1]
controls    <- args[2]
model       <- args[3]
n_cases     <- args[4]
n_vars      <- args[5]
output_file <- args[6] 

# Read and join cases and controls
cases_controls <- dplyr::left_join(
  readr::read_tsv(cases),
  readr::read_tsv(controls)
)

# Function
test_burden <- function(dat, model = 'DOM') {
  # Add a column for the model
  dat <- dplyr::mutate(dat, model = model)

  # Get the 2 x2 cells based on the model
  # Fixed cohort size
  # 0 if count in controls is NA
  dat <- dplyr::mutate(
    dat,
    a = ifelse(model == 'DOM', CASE_COUNT_DOM, CASE_COUNT_REC),
    c = ifelse(model == 'DOM', CONTROL_COUNT_DOM, CONTROL_COUNT_REC),
    c = ifelse(is.na(c), 0, c),
    b = na.omit(unique(dat$CASE_SIZE)),
    d = na.omit(unique(dat$CONTROL_SIZE))
  )

  # Remove rows with 0 in a
  dat <- dplyr::filter(dat, a >= as.integer(n_cases))
  dat <- dplyr::filter(dat, CASE_NVAR >= as.integer(n_vars))
  
  if (nrow(dat) == 0) {
    return(data.frame())
  } else {
    # Perform Fisher's exact test, by gene
    dat <- dplyr::group_split(dat, gene)
    res <- purrr::map_df(
      dat,
      ~{
          # Create a 2 x 2 table
          mat <- with(
            .x,
            rbind(
            c(a, b - a),
            c(c, d - c)
          )
          )

          # clean counts table
          counts <- dplyr::select(
            .x,
            model,
            category,
            gene,
            nvar = CASE_NVAR,
            case_count = a,
            case_size = b,
            control_count = c,
            control_size = d
          )

          # Apply test
          tst <- broom::tidy(
            fisher.test(
              mat
            )
          )
          dplyr::bind_cols(counts, tst)
        }
    )

    # Adjust for multiple testing
    res <- dplyr::mutate(res, p.adj = stats::p.adjust(p.value, method = 'BH'))
  }
  
  return(res)
}

# Call test_burden
gene_burden <- test_burden(
  cases_controls,
  model = model
)

# Write output
readr::write_tsv(
  gene_burden,
  output_file
)
