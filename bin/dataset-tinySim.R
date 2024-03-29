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

# Source functions
suppressMessages({
    source("io.R")
})

#' Install splatter from GitHub using remotes
install_splatter <- function() {
    if (!requireNamespace("splatter", quietly = TRUE)) {
        message("Installing splatter...")
        remotes::install_github(
            "Oshlack/splatter@v1.25.1",
            dependencies = FALSE
        )
    } else {
        message("splatter already installed")
    }
}

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

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    out_file <- args[["--out-file"]]

    install_splatter()

    output <- simulate_dataset()
    print(output)
    write_h5ad(output, out_file, X_name = "counts")
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
