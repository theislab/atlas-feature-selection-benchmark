#!/usr/bin/env Rscript

"
Select features using the singleCellHaystack

Usage:
    method-singleCellHaystack.R --out-file=<path> <file>

Options:
    -h --help               Show this screen.
    --out-file=<path>       Path to output file.
    -n --n-features=<int>   Number of features to select [default: 2000].
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Install singleCellHaystack from CRAN using remotes
install_singleCellHaystack <- function() {
    if (!requireNamespace("singleCellHaystack", quietly = TRUE)) {
        message("Installing singleCellHaystack...")
        remotes::install_cran(
            "singleCellHaystack@0.3.4",
            dependencies = FALSE
        )
    } else {
        message("singleCellHaystack already installed")
    }
}

#' Select features using the singleCellHaystack method
#'
#' @param seurat Seurat object
#' @param n_features Number of features to select
#'
#' @returns data.frame containing the selected features
select_singleCellHaystack_features <- function(seurat, n_features) {

    message("Normalising data...")
    seurat <- Seurat::NormalizeData(seurat)

    message("Selecting highly variable genes...")
    seurat <- Seurat::FindVariableFeatures(seurat)

    message("Scaling data...")
    seurat <- Seurat::ScaleData(seurat)

    message("Calculating PCA...")
    seurat <- Seurat::RunPCA(seurat, features = Seurat::VariableFeatures(seurat))

    message("Selecting ", n_features, " singleCellHaystack features...")
    results <- singleCellHaystack::haystack(seurat, coord = "pca")
    top_results <- singleCellHaystack::show_result_haystack(results, n = n_features)
    top_results$Features <- rownames(top_results)

    return(top_results)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- args[["--n-features"]]

    install_singleCellHaystack()

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
    score <- select_singleCellHaystack_features(seurat, n_features)
    message("Writing output to '", out_file, "'...")
    write.table(
        score,
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
