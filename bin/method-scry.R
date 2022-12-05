#!/usr/bin/env Rscript

"
Select features using the deviance method from scry

Usage:
    method-scry.R --out-file=<path> [options] <file>

Options:
    -h --help                 Show this screen.
    -o --out-file=<path>      Path to output file.
    -n --n-features=<int>     Number of features to select [default: 1000].
" -> doc

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the deviance method from scry
#'
#' @param sce SingleCellExperiment object
#' @param n_features Number of features to select
#'
#' @returns The function output
select_features_scry <- function(sce, n_features) {

    message("Selecting features...")
    sce <- scry::devianceFeatureSelection(sce, nkeep = n_features)
    output <- SummarizedExperiment::rowData(sce)
    output$Feature <- rownames(sce)

    return(output)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- args[["--n-features"]]

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
    print(input)
    output <- select_features_scry(input, n_features)
    message("Writing output to '", out_file, "'...")
    write.table(
        output,
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
