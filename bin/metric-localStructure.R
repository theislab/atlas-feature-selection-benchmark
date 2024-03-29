#!/usr/bin/env Rscript

"
Evaluate integration using the Seurat local structure metric

Usage:
    metric-mixing.R --dataset=<str> --method=<str> --integration=<str> --exprs=<file> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<file>       Path to H5AD file containing the expression matrix.
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

#' Calculate the Seurat local structure metric for an integrated dataset
#'
#' @param seurat Seurat object containing the integrated dataset
#'
#' @returns The local structure metric score
calculate_localStructure <- function(seurat) {
    seurat <- Seurat::NormalizeData(seurat)

    neighbors <- 100
    min_batch_size <- min(table(seurat[["Batch"]]))
    if (neighbors > min_batch_size) {
        warning(
            "some batches have less than 'neighbors' cells, ",
            "setting 'neighbors' to half the smallest batch (",
            floor(min_batch_size / 2), ")"
        )
        neighbors <- floor(min_batch_size / 2)
    }

    message("Calculating local structure scores...")
    cell_scores <- Seurat::LocalStruct(
        seurat,
        "Batch",
        neighbors = neighbors,
        reduction = "emb",
        reduced.dims = seq_len(ncol(seurat[["emb"]])),
        orig.dims = seq_len(ncol(seurat[["emb"]]))
    )

    message("Calculating final score...")
    score <- mean(unlist(cell_scores))

    return(score)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    dataset <- args[["--dataset"]]
    method <- args[["--method"]]
    integration <- args[["--integration"]]
    exprs_file <- args[["--exprs"]]
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
    message("Reading expression data from '", exprs_file, "'...")
    exprs <- read_h5ad(
        exprs_file,
        X_name = "counts",
        uns    = FALSE,
        varm   = FALSE,
        obsm   = FALSE,
        varp   = FALSE,
        obsp   = FALSE
    )
    SingleCellExperiment::reducedDim(exprs, "emb") <- SingleCellExperiment::reducedDim(input, "X_emb")
    message("Converting to Seurat object...")
    seurat <- SeuratObject::as.Seurat(exprs, counts = "counts", data = NULL)
    message("Read data:")
    print(seurat)
    score <- calculate_localStructure(seurat)
    output <- format_metric_results(
        dataset,
        method,
        integration,
        "IntegrationBio",
        "LocalStructure",
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
