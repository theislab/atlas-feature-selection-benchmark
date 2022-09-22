#!/usr/bin/env Rscript

"
Select features using the modelGeneVar method from scran

Usage:
    method-osca.R --out-file=<path> [options] <file>

Options:
    -h --help                 Show this screen.
    -o --out-file=<path>      Path to output file.
    -b --batch                Flag wether to use 'Batch' as block
    -n --n_features=<integer> Number of features to select
" -> doc

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the NBumi method
#'
#' @param input SingleCellExperiment object containing the integrated dataset
#'
#' @returns DataFrame containing the selected features
select_features_scran <- function(input, batch, n_features) {
  
  message("Selecting features using scran...")
  
  sf <- scater::librarySizeFactors(input)
  input <- scater::logNormCounts(input, sf)
  
  if (batch) {
    block = input$Batch
  } else {
    block = NULL
  }
  
  if (is.null(n_features)) {
    n <- NULL
  } else {
    n <- as.integer(n_features)
  }
  
  result <- scran::modelGeneVar(input, block = block)
  
  result <- scran::getTopHVGs(result, n = n)
  
  result <- data.frame(
    "Feature" = result, row.names = result
  )
  
  return(result)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    batch <- args[["--batch"]]
    n_features <- args[["--n_features"]]

    message("Reading data from '", file, "'...")
    input <- read_h5ad(
            file,
            X_name = "counts",
            uns    = FALSE,
            varm   = FALSE,
            obsm   = "X_emb",
            varp   = FALSE,
            obsp   = FALSE
        )
    print(input)
    score <- select_features_scran(input, batch, n_features)
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
