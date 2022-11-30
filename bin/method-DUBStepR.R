#!/usr/bin/env Rscript

"
Select features using Determining the Underlying Basis using Step-wise Regression (DUBStepR) method.
It is designed to select features based on correlations that help the clustering of cell types.

Usage:
    method-DUBStepR.R --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Install DUBStepR from GitHub using remotes
install_DUBStepR <- function() {
    if (!requireNamespace("DUBStepR", quietly = TRUE)) {
        message("Installing DUBStepR...")
        remotes::install_github(
            "prabhakarlab/DUBStepR@76aa39485742d5f5bcfb86346a8a1784ee08f6b9",
            dependencies = FALSE
        )
    } else {
        message("DUBStepR already installed")
    }
}

#' Select features using DUBStepR
#'
#' @param input Seurat object
#'
#' @returns DataFrame containing the selected features.
select_features_dubstepr <- function(seurat) {
    message("Selecting DUBStepR features...")

    message("Normalising data...")
    seurat <- Seurat::NormalizeData(seurat)
    message("Selecting features...")
    results_list <- DUBStepR::DUBStepR(GetAssayData(seurat, slot = "data"))
    result <- results_list$corr.info[results_list$optimal.feature.genes, ]
    colnames(result)[1] <- "Feature"

    return(result)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    install_DUBStepR()

    message("Reading data from '", file, "'...")
    input <- read_h5ad(
        file,
        X_name = "counts",
        uns    = FALSE,
        varm   = FALSE,
        obsm   = FALSE,
        varp   = FALSE,
        obsp   = FALSE
    )

    message("Converting to Seurat object...")
    # Store dummy logcounts for Seurat's conversion function
    SingleCellExperiment::logcounts(input) <- SingleCellExperiment::counts(input)
    seurat <- SeuratObject::as.Seurat(input)
    message("Read data:")
    print(seurat)
    features <- select_features_dubstepr(seurat)
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
