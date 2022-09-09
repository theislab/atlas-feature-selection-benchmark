#!/usr/bin/env Rscript

"
Render an RMarkdown document

Usage:
    render-rmarkdown.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --params=<str>        Comma separated list of key=value parameter pairs
    --out-file=<path>     Path to output file.
" -> doc

#' Render an RMarkdown document
#'
#' @param document Path to the RMarkdown document
#' @param out_file Path to the output file
#' @param params_list list containing parameters passed to the document
#'
#' @returns out_file, invisibly
render_rmarkdown <- function(document, out_file, params_list) {
    message(
        "Rendering '", document,
        "'' to '", out_file,
        "'' with the following parameters:"
    )
    for (item in names(params_list)) {
        message(item, ": ", params_list[[item]])
    }
    rmarkdown::render(
        document,
        output_file = out_file,
        params      = params_list
    )

    invisible(out_file)
}

#' Make a parameters list
#'
#' @param params_str A string representing a set of parameters
#'
#' @returns parameters list
make_params_list <- function(params_str) {
    if (params_str == "") {
        return(list())
    }

    params_split <- stringr::str_split(params_str, ",")[[1]]
    params_split <- stringr::str_split(params_split, "=")

    params_list <- list()
    for (key_value in params_split) {
        params_list[[key_value[1]]] <- key_value[2]
    }

    return(params_list)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    params_str <- args[["--params"]]
    out_file <- args[["--out-file"]]

    params_list <- make_params_list(params_str)
    render_rmarkdown(file, out_file, params_list)
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
