#!/usr/bin/env Rscript

"
Download the scEiaD dataset

Usage:
    dataset-scEiaD.R --out-file=<path> [options]
Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
})

# Source functions
suppressMessages({
    source("io.R")
})

#' Get scEiaD dataset
#'
#' @returns SingleCellExperiment with the scEiaD dataset
get_scEiaD <- function() {
    options(timeout = 3600) # Increase download timeout to 1 hour
    temp_file <- tempfile(fileext = ".Rdata")
    on.exit(file.remove(temp_file))

    url <- "http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/scEiaD_all_seurat_v3.Rdata"
    message("Downloading dataset from ", url, "...")
    result <- download.file(url, temp_file)
    if (result != 0) {
        stop("Failed to download dataset")
    }

    message("Loading dataset...")
    load(temp_file)

    message("Creating SingleCellExperiment object...")
    sce <- SingleCellExperiment(
        assays = list(counts = GetAssayData(scEiaD, slot = "counts")),
        colData = scEiaD[[]]
    )
    rm(scEiaD)
    gc()

    message("Selecting Human cells...")
    is_human <- sce$organism == "Homo sapiens"

    message("Selecting Eye cells...")
    is_eye <- sce$Organ == "Eye"

    message("Selecting Tissue cells...")
    is_tissue <- sce$Source == "Tissue"

    message("Selecting cells with labels...")
    has_label <- !is.na(sce$CellType_predict)

    message("Selecting non-doublet cells...")
    non_doublet <- sce$Doublet == "FALSE"

    message("Subsetting to selected cells...")
    selected <- is_human & is_eye & is_tissue & has_label & non_doublet
    sce <- sce[, selected]

    return(sce)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    out_file <- args[["--out-file"]]

    output <- get_scEiaD()
    print(output)
    write_h5ad(output, out_file, X_name = "counts")
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
