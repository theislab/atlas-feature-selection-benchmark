#!/usr/bin/env Rscript

"
Evaluate integration using the Cell-Specific Mixing Score (CMS) metric

Usage:
    metric-cms.R --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
" -> doc

# Source functions
suppressMessages({
    source("io.R")
    source("metrics.R")
})

#' Calculate the Cell-Specific Mixing Score metric for an integrated dataset
#'
#' @param sce SingleCellExperiment object containing the integrated dataset
#'
#' @returns The CMS metric score
calculate_CMS <- function(sce) {
    set.seed(1)

    message("Calculating cell mixing scores...")
    n_cores <- parallelly::availableCores()
    message("Using ", n_cores, " cores")
    sce <- CellMixS::cms(
        sce,
        k       = 200,
        group   = "Batch",
        dim_red = "X_emb",
        n_dim   = ncol(SingleCellExperiment::reducedDim(sce, "X_emb")),
        BPPARAM = BiocParallel::MulticoreParam(workers = n_cores)
    )

    message("Calculating final score...")
    # Score is 1 minus the proportion of cells with a CMS p-value < 0.1
    cms_scores <- SummarizedExperiment::colData(sce)$cms
    score <- 1 - mean(cms_scores < 0.1)

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
    message("Read data:")
    print(input)
    score <- calculate_CMS(input)
    output <- format_metric_results(
        dataset,
        method,
        integration,
        "IntegrationBatch",
        "CMS",
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
