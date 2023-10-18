#!/usr/bin/env Rscript

"
Download the Human Transcription Factors dataset

Usage:
    dataset-human-tfs.R --out-file=<path> [options]

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

#' Get the Human Transcription Factors gene set
#'
#' Download the genes from the website, select TFs and use {biomaRt} to add
#' mouse IDs/symbols
#'
#' @returns data.frame with Human transcription factor genes
get_human_tfs <- function() {
    url <- "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"

    message("Reading database file from '", url, "'...")
    database <- read.csv(url, row.names = 1)

    message("Selecting TFs...")
    is_tf <- database$Is.TF. == "Yes"
    human_tfs <- database[is_tf, c("Ensembl.ID", "HGNC.symbol")]
    colnames(human_tfs) <- c("ENSEMBL", "Gene")
    human_tfs$Species <- "Human"

    message("Getting human mart...")
    human_mart <- biomaRt::useEnsembl(
        "ensembl", "hsapiens_gene_ensembl",
        version = "105"
    )

    message("Getting mouse mart...")
    mouse_mart <- biomaRt::useEnsembl(
        "ensembl", "mmusculus_gene_ensembl",
        version = "105"
    )

    message("Getting mouse genes...")
    mapping <- biomaRt::getLDS(
        attributes  = "ensembl_gene_id",
        filters     = "ensembl_gene_id",
        values      = human_tfs$ENSEMBL,
        mart        = human_mart,
        attributesL = c("ensembl_gene_id", "mgi_symbol"),
        martL       = mouse_mart
    )
    mouse_tfs <- data.frame(
        ENSEMBL = mapping$Gene.stable.ID.1,
        Gene    = mapping$MGI.symbol,
        Species = "Mouse"
    )

    message("Creating final data.frame...")
    tfs <- rbind(human_tfs, mouse_tfs)
    tfs <- tfs[!duplicated(tfs), ]

    return(tfs)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    out_file <- args[["--out-file"]]

    output <- get_human_tfs()
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
