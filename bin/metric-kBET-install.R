#!/usr/bin/env Rscript

"
Install kBET from GitHub

Usage:
    metric-kBET-install.R

Options:
    -h --help             Show this screen.
" -> doc

#' Install kBET from GitHub using remotes
install_kBET <- function() {
    if (!requireNamespace("kBET", quietly = TRUE)) {
        message("Installing kBET...")
        remotes::install_github(
            "theislab/kBET@a10ffeaa31da83e4305dfe85cd0adfcebeee721e",
            dependencies = FALSE
        )
    } else {
        message("kBET already installed")
    }
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)
    install_kBET()
}

if (sys.nframe() == 0) {
    main()
}
