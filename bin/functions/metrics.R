#' Format metric results
#'
#' @param dataset The name of the dataset the metric has been calculated for
#' @param method The name of the method the metric has been calculated for
#' @param integration The name of the integration the metric has been calculated
#' for
#' @param metric_type The type of the metric, either 'Integration',
#' 'Classification', 'Mapping' or 'Unseen'
#' @param metric The name of the metric that has been calculated
#' @param value The value of the calculated metric
#'
#' @return data.frame containing the formatted results
format_metric_results <- function(dataset, method, integration, metric_type, metric, value) {
    if (!(metric_type %in% c("Integration", "Classification", "Mapping", "Unseen"))) {
        stop("'metric_type' must be one of 'Integration', 'Classification', 'Mapping' or 'Unseen'")
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
