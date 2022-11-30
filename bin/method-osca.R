#!/usr/bin/env Rscript

"
Select features using method described in 'Orchestrating Single-Cell Analysis with Bioconductor'

Usage:
    method-osca.R --out-file=<path> [options] <file>

Options:
    -h --help                 Show this screen.
    -o --out-file=<path>      Path to output file.
    -n --n_features=<int>     Number of features to select [default: 1000].
" -> doc

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the OSCA method
#'
#' @param input SingleCellExperiment object
#' @param n_features The number of features to select
#'
#' @returns DataFrame containing the selected features
select_features_osca <- function(input, n_features) {

  message("Selecting features using the OSA method...")

  message("Performing batch-aware normalisation...")
  input <- batchelor::multiBatchNorm(input, batch = input$Batch)

  message("Moodelling features variance in a batch-aware way...")
  feature_stats <- scran::modelGeneVar(input, block = input$Batch)

  message("Selecting top variance features...")
  top_hvgs <- scran::getTopHVGs(feature_stats, n = n_features)

  result <- feature_stats[top_hvgs, ]
  result$Feature <- top_hvgs
  result$per.block <- NULL

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
    output <- select_features_osca(input, n_features)
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
