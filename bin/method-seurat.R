#!/usr/bin/env Rscript

"
Select features using seurat

Usage:
    method-seurat.R --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    -o --out-file=<path>    Path to output file.
    -n --n-features=<int>   Number of features to select [default: 1000].
    -f --flavor=<str>       Flavor of feature selection method. One of: seurat, cell_ranger, seurat_v3 [default: 'vst'].
" -> doc


# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using seurat
#'
#' @param seurat Seurat object
#'
#' @returns DataFrame containing the selected features
select_seurat_features <- function(seurat, n_features, flavor) {

    message("Selecting seurat features...")

    result <- FindVariableFeatures(
        seurat,
        selection.method = flavor,
        nfeatures = n_features
    )
    selected_features <- data.frame(Feature= VariableFeatures(result))

    return(selected_features)
}

#' The main script function
main <- function() {

    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- args[["--n-features"]]
    flavor <- args[['--flavor']]

    message("Reading data from '", file, "'...")
    input <- read_h5ad(
            file,
            X_name = 'counts',
            uns    = FALSE,
            varm   = FALSE,
            obsm   = "X_emb",
            varp   = FALSE,
            obsp   = FALSE
        )

    message("Converting to Seurat object...")
    # Store dummy logcounts for Seurat's conversion function
    SingleCellExperiment::logcounts(input) <- SingleCellExperiment::counts(input)
    seurat <- SeuratObject::as.Seurat(input)
    message("Read data:")
    print(seurat)
    score <- select_seurat_features(seurat, n_features, flavor)
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
