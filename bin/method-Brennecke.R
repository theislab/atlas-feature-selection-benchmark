#!/usr/bin/env Rscript

"
Select features using the excess CV method from Brennecke et al.

Usage:
    method-Brennecke.R --out-file=<path> [options] <file>

Options:
    -h --help                 Show this screen.
    -o --out-file=<path>      Path to output file.
    -n --n_features=<int>     Number of features to select [default: 1000].
" -> doc

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the Brennecke method
#'
#' @param input SingleCellExperiment object
#' @param n_features The number of features to select
#'
#' @returns DataFrame containing the selected features
select_features_Brennecke <- function(input, n_features) {
    message("Selecting features using the OSA method...")

    message("Normalising expression...")
    input <- scuttle::logNormCounts(input)

    message("Modelling CV...")
    feature_stats <- scran::modelGeneCV2(input)

    message("Selecting top CV features...")
    top_hvgs <- scran::getTopHVGs(feature_stats, var.field = "ratio", n = n_features)

    result <- feature_stats[top_hvgs, ]
    result$Feature <- top_hvgs

    return(result)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- args[["--n_features"]]

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
    output <- select_features_Brennecke(input, n_features)
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
