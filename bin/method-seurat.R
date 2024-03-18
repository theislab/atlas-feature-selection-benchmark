#!/usr/bin/env Rscript

"
Select features using seurat

Usage:
    method-seurat.R --out-file=<path> [options] <file>

Options:
    -h --help               Show this screen.
    -o --out-file=<path>    Path to output file.
    -n --n-features=<int>   Number of features to select [default: 2000].
    -m --method=<str>       Feature selection method to use. One of: vst, mvp, disp, sct [default: 'vst'].
" -> doc


# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
})

# Source functions
suppressMessages({
    source("io.R")
})

#' Select features using Seurat
#'
#' @param seurat Seurat object
#' @param n_features The number of features to select
#' @param method The feature selection method to use
#'
#' @returns DataFrame containing the selected features
select_seurat_features <- function(seurat, n_features,
    method = c("vst", "mvp", "disp", "sct")) {
    method <- match.arg(method)

    method <- switch(method,
        vst  = "vst",
        mvp  = "mean.var.plot",
        disp = "dispersion",
        sct  = "sctransform"
    )

    message(
        "Selecting ", n_features,
        " features using the Seurat '",
        method, "' method..."
    )

    if (method == "sctransform") {
        seurat <- SCTransform(
            seurat,
            assay = DefaultAssay(seurat),
            variable.features.n = n_features,
            do.correct.umi = FALSE,
            conserve.memory = TRUE,
            verbose = TRUE
        )
    } else {
        seurat <- NormalizeData(seurat)
        seurat <- FindVariableFeatures(
            seurat,
            selection.method = method,
            nfeatures = n_features
        )
    }

    selected_features <- data.frame(Feature = VariableFeatures(seurat))

    return(selected_features)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]
    n_features <- args[["--n-features"]]
    method <- args[["--method"]]

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

    message("Converting to Seurat object...")
    # Store dummy logcounts for Seurat's conversion function
    seurat <- SeuratObject::as.Seurat(input, counts = "counts", data = NULL)
    message("Read data:")
    print(seurat)
    output <- select_seurat_features(seurat, n_features, method)
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
