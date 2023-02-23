#!/usr/bin/env Rscript

"
Evaluate cell label classification using the accuracy metric for the rarest cell label

Usage:
    metric-rareAccuracy.R --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
" -> doc

# Source functions
suppressMessages({
    source("functions.R")
})

#' Calculate classification accuracy for the rarest cell label
#'
#' @param labels data.frame containing real and predicted labels
#'
#' @returns Accuracy for the rarest cell type
calculate_rare_accuracy <- function(labels) {
    message("Identify rare label...")
    label_freqs <- table(labels$Label)
    rare_label <- names(label_freqs)[which.min(label_freqs)]

    message(
        "Calculating classification accuracy for '", rare_label, "' label..."
    )
    is_rare <- labels$Label == rare_label
    n_correct <- sum(labels$PredLabel[is_rare] == rare_label)
    score <- n_correct / sum(is_rare)

    return(score)
}

#' The main script function
main <- function() {
    args <- docopt::docopt(doc)

    file <- args[["<file>"]]
    dataset <- args[["--dataset"]]
    method <- args[["--method"]]
    integration <- args[["--integration"]]
    out_file <- args[["--out-file"]]

    message("Reading data from '", file, "'...")
    input <- readr::read_tsv(
        file,
        col_types = readr::cols(
            .default  = readr::col_double(),
            ID        = readr::col_character(),
            Unseen    = readr::col_logical(),
            Label     = readr::col_character(),
            PredLabel = readr::col_character()
        )
    )
    message("Removing unseen populations...")
    input <- input[!input$Unseen, ]
    message("Read data:")
    print(input)
    score <- calculate_rare_accuracy(input)
    output <- format_metric_results(
        dataset,
        method,
        integration,
        "Classification",
        "RareAccuracy",
        score
    )
    print(output)
    message("Writing output to '", out_file, "'...")
    readr::write_tsv(output, out_file)
    message("Done!")
}

if (sys.nframe() == 0) {
    main()
}
