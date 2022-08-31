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

#' Read a H5AD file to a SingleCellExperiment object
#'
#' @param file Path to the H5AD file
#' @param ... Additional parameters passed to `zellkonverter::AnnData2SCE`
#'
#' @return `SingleCellExperiment` object
read_h5ad <- function(file, ...) {
    anndata <- reticulate::import("anndata")

    message("Reading AnnData from ", file, "...")
    adata <- anndata$read_h5ad(file)

    message("Converting to SingleCellExperiment...")
    sce <- zellkonverter::AnnData2SCE(
        adata,
        verbose = TRUE,
        ...
    )

    return(sce)
}

#' Format metric results
#'
#' @param dataset The name of the dataset the metric has been calculated for
#' @param method The name of the method the metric has been calculated for
#' @param integration The name of the integration the metric has been calculated
#' for
#' @param metric_type The type of the metric, either "Integration" or
#' "Classification"
#' @param metric The name of the metric that has been calculated
#' @param value The value of the calculated metric
#'
#' @return data.frame containing the formatted results
format_metric_results <- function(dataset, method, integration, metric_type, metric, value) {
    if (!(metric_type %in% c("Integration", "Classification"))) {
        stop("'metric_type' must be one of 'Integration' or 'Classification'")
    }

    if (value < 0 || value > 1) {
        stop("'score' must be between 0 and 1")
    }

    message("Formatting metric results...")
    data.frame(
        "Dataset"     = dataset,
        "Method"      = method,
        "Integration" = integration,
        "Type"        = metric_type,
        "Metric"      = metric,
        "Value"       = value
    )
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
    # Store dummy logcounts for Seurat's conversion function
    SingleCellExperiment::logcounts(input) <- SingleCellExperiment::counts(input)
    SingleCellExperiment::reducedDimNames(input) <- "emb"
    seurat <- SeuratObject::as.Seurat(input)
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
