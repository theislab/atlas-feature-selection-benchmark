#!/usr/bin/env Rscript

"
Simulate a tiny dataset

Usage:
    dataset-tinySim.R --out-file=<path> [options]
Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

#' Simulate a tiny dataset
#'
#' @returns SingleCellExperiment with simulated data
simulate_dataset <- function() {
    message("Simulating data...")
    sim <- splatter::splatSimulateGroups(
        nGenes       = 10000,
        batchCells   = c(100, 100, 100, 100),
        batch.facLoc = c(0.15, 0.15, 0.15, 0.20),
        group.prob   = c(0.30, 0.25, 0.20, 0.15, 0.10),
        de.prob      = c(0.05, 0.10, 0.20, 0.30, 0.10),
        lib.loc      = 10,
        seed         = 1,
        verbose      = TRUE
    )

    message("Performing quick cell filtering...")
    sim <- scuttle::quickPerCellQC(sim)

    return(sim)
}

#' Write a SingleCellExperiment to a H5AD file
#'
#' @param sce The SingleCellExperiment object to write
#' @param file The file to write to
#'
#' @return `file` invisibly
write_h5ad <- function(sce, file) {
    anndata <- reticulate::import("anndata")

    message("Converting to AnnData...")
    adata <- zellkonverter::SCE2AnnData(
        sce,
        X_name  = "counts",
        verbose = TRUE
    )

    message("Writing AnnData to '", file, "'...")
    adata$write_h5ad(file)

    invisible(file)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    out_file <- args[["--out-file"]]

    output <- simulate_dataset()
    print(output)
    write_h5ad(output, out_file)
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
