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
