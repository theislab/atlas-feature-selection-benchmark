#!/usr/bin/env Rscript

"
Combine and summarise metric scores

Usage:
    combine-metrics.R --out-file=<path> [options] <file>...

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

#' Combine metrics output files
#'
#' @param files Vector of metric output file paths
#'
#' @returns Combined metrics tibble
combine_metrics <- function(files) {
    message("Reading metrics from ", length(files), " files...")
    metrics <- readr::read_tsv(
        files,
        col_types = readr::cols(
            .default = readr::col_character(),
            Value    = readr::col_double()
        )
    )

    message("Scaling metric scores...")
    method_names <- unique(metrics$Method)
    n_random <- sum(stringr::str_detect(method_names, "random"))
    mode <- "range"
    if (n_random <= 1) {
        warning(
            "Scores for less than two random gene sets found. ",
            "All values will be used for scaling."
        )
        mode <- "all"
    }
    metrics <- metrics |>
        dplyr::group_by(Dataset, Metric, Integration) |>
        dplyr::mutate(
            ScaledScore = scale_values(
                Value,
                stringr::str_detect(Method, "random"),
                mode = mode
            )
        )

    return(metrics)
}

#' Scale values
#'
#' Perform z-score scaling on a set of values
#'
#' @param values Vector of values to scale
#' @param is_range Logical vector setting which values to use to define the
#' reference range (if `mode == "range"`)
#' @param mode The scaling mode to use, either "range" to only use the values
#' specified by `is_range` or "all" to use all values
#'
#' @details
#' If `mode == "range"` only a subset of values corresponding to a reference
#' range are used to calculate the mean and standard deviation which is then
#' used to scale all values. For `mode == "all"` normal scaling using all values
#' is performed.
#'
#' If the standard deviation is zero (all scores are the same) then a score of
#' zero is returned with a warning (rather than `NaN`)
#'
#' @returns Vector of scaled values
scale_values <- function(values, is_range, mode = c("range", "all")) {
    mode <- match.arg(mode)

    mu <- switch(mode,
        range = mean(values[is_range]),
        all   = mean(values)
    )

    sigma <- switch(mode,
        range = sd(values[is_range]),
        all   = sd(values)
    )

    if (sigma == 0) {
        warning("No variation in scores, setting scaled scores to 0")
        return(rep(0, length(values)))
    }

    scale(values, center = mu, scale = sigma)[, 1]
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    files <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    message("Reading data from '", file, "'...")
    input <- read_data(file)
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
