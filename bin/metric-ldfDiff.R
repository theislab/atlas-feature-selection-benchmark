#!/usr/bin/env Rscript

"
Evaluate integration using the Difference in Local Density Factor (ldfDiff) metric

Usage:
    metric-ldfDiff.R --dataset=<str> --method=<str> --integration=<str> --exprs=<file> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<file>       Path to H5AD file containing the expression matrix.
    --out-file=<path>    Path to output file.
" -> doc

# Source functions
suppressMessages({
    source("io.R")
    source("metrics.R")
})

#' Calculate the Difference in Local Density Factor metric for an integrated dataset
#'
#' @param input SingleCellExperiment object containing the integrated dataset
#' @param exprs SingleCellExperiment object containing the expression matrix
#'
#' @returns The ldfDiff metric score
calculate_ldfDiff <- function(input, exprs) {
    set.seed(1)

    n_dim <- ncol(SingleCellExperiment::reducedDim(input, "X_emb"))

    message("Calculating batch PCAs...")
    batches <- sort(unique(exprs[["Batch"]]))
    batch_objects <- lapply(batches, function(.batch) {
        message("Batch '", .batch, "'...")
        batch_sce <- exprs[, exprs[["Batch"]] == .batch]
        batch_sce <- scuttle::logNormCounts(batch_sce)
        batch_sce <- scater::runPCA(
            batch_sce,
            exprs_values = "logcounts",
            ncomponents  = 10,
            name         = "X_pca" # Use a different name just to make sure it isn't recalculated
        )
    })
    names(batch_objects) <- batches

    # Use k=75 unless one of the batches is smaller than that
    # (should only happen for the test dataset)
    k <- min(sapply(batch_objects, ncol)) - 1
    if (k < 75) {
        warning(
            "k was set to ", k,
            " because one of the batches has fewer than 76 cells"
        )
    } else {
        k <- 75
    }
    message("Calculating cell ldfDiff scores...")
    input <- CellMixS::ldfDiff(
        sce_pre_list = batch_objects,
        sce_combined = input,
        group        = "Batch",
        k            = k,
        dim_red      = "X_pca",
        dim_combined = "X_emb",
        n_dim        = n_dim
    )

    message("Calculating final ldfDiff score...")
    # Scores closer to 0 are better so use absolute values
    scores <- abs(SummarizedExperiment::colData(input)$diff_ldf)
    if (any(is.na(scores))) {
        message(
            "Warning: Ignoring ", sum(is.na(scores)), " cells with NA ldfDiff scores"
        )
        scores <- scores[!is.na(scores)]
    }
    # Score are unbounded so we set any scores greater than 1 to 1
    if (any(scores > 1)) {
        warning(
            sum(scores > 1), " absolute cell scores (",
            round(sum(scores > 1) / ncol(input) * 100, 2),
            "%) are greater than 1. Setting these to 1."
        )
        scores[scores > 1] <- 1
    }
    # Take the mean of the scores and subtract it from 1 so higher is better
    score <- 1 - mean(scores)

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
    message("Read data:")
    print(input)
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
    print(exprs)
    score <- calculate_ldfDiff(input, exprs)
    output <- format_metric_results(
        dataset,
        method,
        integration,
        "IntegrationBio",
        "ldfDiff",
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
