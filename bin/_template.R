#!/usr/bin/env Rscript

"
An R script template

Usage:
    template.R --out-file=<path> [options] <file>
Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

# Load libraries
suppressPackageStartupMessages({

})

#' A function that performs analysis
#'
#' @param input The function input
#'
#' @returns The function output
run_analysis <- function(input) {

    output = input

    return(input)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file      <- args[["<file>"]]
    out_file  <- args[["--out-file"]]

    message("Reading data from '", file, "'...")
    # input <- read_data(file)
    message("Read data:")
    # print(input)
    output <- run_analysis(input)
    message("Writing output to '", out_file, "'...")
    # write_output(output, out_file)
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
