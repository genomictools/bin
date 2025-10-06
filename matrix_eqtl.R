#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

snps_file  <- args[1]
traits_file<- args[2]
cvrt_file  <- args[3]
output     <- args[4]

# modelANOVA or modelLINEAR or modelLINEAR_CROSS
useModel <- MatrixEQTL::modelLINEAR

snps <- MatrixEQTL::SlicedData$new()
snps$LoadFile( snps_file )

traits <- MatrixEQTL::SlicedData$new()
traits$LoadFile( traits_file )

cvrt <- MatrixEQTL:: SlicedData$new()
cvrt$LoadFile( cvrt_file )

# snps$fileDelimiter = "\t"      # the TAB character
# snps$fileOmitCharacters = "NA" # denote missing values
# snps$fileSkipRows = 1          # one row of column labels
# snps$fileSkipColumns = 1       # one column of row labels
# snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows

# Call Matrix eQTL
me <- MatrixEQTL::Matrix_eQTL_engine(
    snps = snps,
    gene = traits,
    cvrt = cvrt,
    output_file_name = output,
    pvOutputThreshold = 1e-2,
    errorCovariance = numeric(),
    useModel = useModel,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)
