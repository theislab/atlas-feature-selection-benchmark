#!/usr/bin/env Rscript

"
Evaluate integration using the Seurat mixing metric

Usage:
    metric-mixing.R --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
})

# Source functions
suppressMessages({
    source("io.R")
    source("metrics.R")
})

#' Calculate the Seurat mixing metric for an integrated dataset
#'
#' @param seurat Seurat object containing the integrated dataset
#'
#' @returns The mixing metric score
calculate_mixing <- function(seurat) {
    max_k <- 300
    if (max_k > ncol(seurat)) {
        warning(
            "'max_k' greater than the number of cells, setting 'max_k' to ",
            ncol(seurat)
        )
        max_k <- ncol(seurat)
    }

    message("Calculating cell mixing scores...")
    cell_scores <- Seurat::MixingMetric(
        seurat,
        "Batch",
        reduction = "emb",
        dims      = ncol(seurat[["emb"]]),
        max.k     = max_k
    )

    message("Calculating final score...")
    # Calculate the mean normalised by the max score as the overall score
    # Subtract from 1 so higher scores are better
    score <- 1 - mean(cell_scores / max_k)

    return(score)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    dataset <- args[["--dataset"]]
    method <- args[["--method"]]
    integration <- args[["--integration"]]
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
    message("Converting to Seurat object...")
    SingleCellExperiment::reducedDimNames(input) <- "emb"
    seurat <- SeuratObject::as.Seurat(input, counts = "counts", data = NULL)
    message("Read data:")
    print(seurat)
    score <- calculate_mixing(seurat)
    output <- format_metric_results(
        dataset,
        method,
        integration,
        "Integration",
        "Mixing",
        score
    )
    print(output)
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
