#!/usr/bin/env Rscript

"
Select features using the scPNMF method

Usage:
    method-scPNMF.R --out-file=<path> [options] <file>

Options:
    -h --help               Show this screen.
    --out-file=<path>       Path to output file.
    -n --n-features=<int>   Number of features to select [default: 200].
" -> doc

# Source functions
suppressMessages({
    source("io.R")
})

#' Install scPNMF from GitHub using remotes
install_scPNMF <- function() {
    if (!requireNamespace("akmedoids", quietly = TRUE)) {
        message("Installing akmedoids...")
        remotes::install_version(
            "akmedoids",
            version = "1.3.0",
            repos = "https://cloud.r-project.org",
            dependencies = FALSE
        )
    } else {
        message("akmedoids already installed")
    }

    if (!requireNamespace("scPNMF", quietly = TRUE)) {
        message("Installing scPNMF...")
        remotes::install_github(
            "JSB-UCLA/scPNMF@47d5b10cb09450255aea9b53ace555a95ab69502",
            dependencies = FALSE
        )
    } else {
        message("scPNMF already installed")
    }
}

#' Select features using the scPNMF method
#'
#' @param sce SingleCellExperiment object
#' @param n_features Number of features to select
#'
#' @returns data.frame containing the selected features
select_scPNMF <- function(sce, n_features) {
    message("Normalising counts...")
    sce <- scuttle::logNormCounts(sce)

    message("Fitting PNMF model...")
    pnmf <- scPNMF::PNMFfun(
        SingleCellExperiment::logcounts(sce),
        K        = 20, # Recommended by vignette
        method   = "EucDist",
        tol      = 1e-4,
        maxIter  = 1000,
        verboseN = TRUE
    )

    message("Selecting bases...")
    W_select <- scPNMF::basisSelect(
        W          = pnmf$Weight,
        S          = pnmf$Score,
        X          = SingleCellExperiment::logcounts(sce),
        toTest     = TRUE,
        toAnnotate = FALSE,
        mc.cores   = 1
    )

    message("Selecting informative genes...")
    genes <- scPNMF::getInfoGene(
        W_select,
        M            = n_features,
        by_basis     = FALSE,
        return_trunW = TRUE,
        dim_use      = NULL
    )

    output <- data.frame(Feature = genes$InfoGene)

    return(output)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- as.numeric(args[["--n-features"]])

    install_scPNMF()

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
    message("Read data:")
    print(input)
    output <- select_scPNMF(input, n_features)
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
