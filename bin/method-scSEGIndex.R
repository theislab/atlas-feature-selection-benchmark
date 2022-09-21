#!/usr/bin/env Rscript

"
Select features using the scSEGIndex method from scMerge

Usage:
    method-scSEGIndex.R --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
" -> doc


# Load libraries
suppressPackageStartupMessages({
    library(scMerge)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the scSEGIndex method
#'
#' @param seurat Seurat object containing the integrated dataset
#'
#' @returns DataFrame containing the selected features
select_scsegindex_features <- function(input) {

    message("Selecting scSEGIndex features...")

    exprs_mat = SummarizedExperiment::assay(input)
    result = scSEGIndex(exprs_mat = exprs_mat)
    result$Feature = rownames(result)

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
    score <- select_scsegindex_features(input)
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
