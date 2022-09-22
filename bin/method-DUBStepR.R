#!/usr/bin/env Rscript

"
Select features using Determining the Underlying Basis using Step-wise Regression (DUBStepR) method.
It is aimed to select features based on correlations that help the clustering of cell types.

Usage:
    method-DUBStepR.R --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
" -> doc


# Load libraries
suppressPackageStartupMessages({
    library(DUBStepR)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using DUBStepR.
#'
#' @param input SingleCellExperiment object containing the integrated dataset.
#'
#' @returns DataFrame containing the selected features.
select_features_dubstepr <- function(input) {

    message("Selecting DUBStepR features...")

    exprs_mat <- SummarizedExperiment::assay(input)
    result0 <- DUBStepR::DUBStepR(input.data = exprs_mat, min.cells = 0.05*ncol(exprs_mat))
    result <- result0$corr.info[result0$optimal.feature.genes, ]
    colnames(result)[1] <- "Feature"

    return(result)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    
    message("Reading data from '", file, "'...")
    input <- read_h5ad(
        file,
        X_name = NULL,
        uns    = FALSE,
        varm   = FALSE,
        obsm   = "X_emb",
        varp   = FALSE,
        obsp   = FALSE
    )
    print(input)
    features <- select_features_dubstepr(input)
    message("Writing output to '", out_file, "'...")
    write.table(
        features,
        file      = out_file,
        quote     = FALSE,
        sep       = "\t",
        row.names = FALSE
    )
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
