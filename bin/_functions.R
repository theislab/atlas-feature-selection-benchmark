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

#' Write a SingleCellExperiment to a H5AD file
#'
#' @param sce The SingleCellExperiment object to write
#' @param file The file to write to#'
#' @param ... Additional parameters passed to `zellkonverter::AnnData2SCE`
#'
#' @return `file` invisibly
write_h5ad <- function(sce, file, ...) {
    anndata <- reticulate::import("anndata")

    message("Converting to AnnData...")
    adata <- zellkonverter::SCE2AnnData(
        sce,
        verbose = TRUE,
        ...
    )

    message("Writing AnnData to '", file, "'...")
    adata$write_h5ad(file)

    invisible(file)
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
    if (!(metric_type %in% c("Integration", "Classification", "Mapping"))) {
        stop("'metric_type' must be one of 'Integration', 'Classification' or 'Mapping'")
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
