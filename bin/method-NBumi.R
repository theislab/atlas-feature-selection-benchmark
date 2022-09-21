#!/usr/bin/env Rscript

"
Select features using the NBumi method from M3Drop

Usage:
    method-NBumi.R --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
" -> doc


# Load libraries
suppressPackageStartupMessages({
    library(M3Drop)
})

# Source functions
suppressMessages({
    source("_functions.R")
})

#' Select features using the NBumi method
#'
#' @param input SingleCellExperiment object containing the integrated dataset
#'
#' @returns DataFrame containing the selected features
select_nbumi_features <- function(input) {

    message("Selecting NBumi features...")

    count_mat <- NBumiConvertData(SingleCellExperiment::counts(input), is.counts=TRUE)
    DANB_fit <- NBumiFitModel(count_mat)
    result <- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.01, suppress.plot=TRUE)
    result$Feature <- rownames(result)

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
            X_name = "counts",
            uns    = FALSE,
            varm   = FALSE,
            obsm   = "X_emb",
            varp   = FALSE,
            obsp   = FALSE
        )
    print(input)
    score <- select_nbumi_features(input)
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
