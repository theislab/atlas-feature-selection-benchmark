#!/usr/bin/env Rscript

"
Download the Tirosh cell cycle genes

Usage:
    dataset-tirosh-genes.R --out-file=<path> [options]

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

#' Get the Tirosh cell cycle genes
#'
#' Download the genes from the scIB repository and use {biomaRt} to add
#' ENSEMBL gene symbols.
#'
#' @returns data.frame with Tirosh cell cycle genes
get_tirosh_genes <- function() {

    # Base URL for the scIB repository
    base_url <- "https://raw.githubusercontent.com/theislab/scib/cce7aaa65ccc18649e141e30464361d1fc67ea77/scib/resources/"

    # Files containing the Tirosh cell cycle genes
    files <- c(
        mouse_G2M = "g2m_genes_tirosh.txt",
        human_G2M = "g2m_genes_tirosh_hm.txt",
        mouse_S   = "s_genes_tirosh.txt",
        human_S   = "s_genes_tirosh_hm.txt"
    )

    message("Reading gene set files...")
    gene_sets <- lapply(files, function(.file) {
        url <- paste0(base_url, .file)
        message("Reading '", url, "'...")
        read.table(url, header = FALSE, col.names = "Gene")
    })
    names(gene_sets) <- names(files)

    message("Getting human gene annotations...")
    # Get the human mart
    human_mart <- biomaRt::useEnsembl("ensembl","hsapiens_gene_ensembl")

    # Get the annotations for the human genes
    human_annot <- biomaRt::getBM(
        attributes = c("ensembl_gene_id","hgnc_symbol"),
        filters    = "hgnc_symbol",
        values     = c(gene_sets$human_G2M$Gene, gene_sets$human_S$Gene),
        mart       = human_mart
    )
    # Remove any duplicates
    human_annot <- human_annot[!duplicated(human_annot), ]

    # Create a vector mapping the gene names to ENSEMBL gene symbols
    human_map <- human_annot$ensembl_gene_id
    names(human_map) <- human_annot$hgnc_symbol

    message("Adding human gene annotations...")
    # Alse set phase and species
    gene_sets$human_G2M$ENSEMBL <- human_map[gene_sets$human_G2M$Gene]
    gene_sets$human_G2M$Phase   <- "G2M"
    gene_sets$human_G2M$Species <- "Human"
    gene_sets$human_S$ENSEMBL   <- human_map[gene_sets$human_S$Gene]
    gene_sets$human_S$Phase     <- "S"
    gene_sets$human_S$Species   <- "Human"

    message("Getting mouse gene annotations...")
    # Get the human mart
    mouse_mart <- biomaRt::useEnsembl("ensembl","mmusculus_gene_ensembl")

    # Get the annotations for the human genes
    mouse_annot <- biomaRt::getBM(
        attributes = c("ensembl_gene_id","mgi_symbol"),
        filters    = "mgi_symbol",
        values     = c(gene_sets$mouse_G2M$Gene, gene_sets$mouse_S$Gene),
        mart       = mouse_mart
    )
    # Remove any duplicates
    mouse_annot <- mouse_annot[!duplicated(mouse_annot), ]

    # Create a vector mapping the gene names to ENSEMBL gene symbols
    mouse_map <- mouse_annot$ensembl_gene_id
    names(mouse_map) <- mouse_annot$mgi_symbol

    message("Adding mouse gene annotations...")
    # Alse set phase and species
    gene_sets$mouse_G2M$ENSEMBL <- mouse_map[gene_sets$mouse_G2M$Gene]
    gene_sets$mouse_G2M$Phase   <- "G2M"
    gene_sets$mouse_G2M$Species <- "Mouse"
    gene_sets$mouse_S$ENSEMBL   <- mouse_map[gene_sets$mouse_S$Gene]
    gene_sets$mouse_S$Phase     <- "S"
    gene_sets$mouse_S$Species   <- "Mouse"

    message("Creating final data.frame...")
    tirosh_genes <- do.call(rbind, gene_sets)

    return(tirosh_genes)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    output <- get_tirosh_genes()
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
