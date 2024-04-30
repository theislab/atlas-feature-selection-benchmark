#' Summarise metrics
#'
#' Calculate type average and overall metric values
#'
#' @param metrics data.frame containing metrics
#' @param baseline_ranges date.frame containing baseline ranges
#'
#' @return data.frame containing summarised metrics
summarise_metrics <- function(metrics, baseline_ranges) {
    metrics |>
        scale_metrics(baseline_ranges) |>
        dplyr::group_by(
            .data$Dataset, .data$Method, .data$Integration, .data$Type
        ) |>
        dplyr::summarise(
            TypeMean = mean(.data$ScaledValue),
            .groups = "drop"
        ) |>
        tidyr::pivot_wider(
            names_from = "Type", values_from = "TypeMean"
        ) |>
        dplyr::mutate(
            Overall = 0.5 * (
                ((1 / 2) * .data$IntegrationBatch) +
                    ((1 / 2) * .data$IntegrationBio)
            ) +
                0.5 * (
                    ((1 / 3) * .data$Mapping) +
                        ((1 / 3) * .data$Classification) +
                        ((1 / 3) * .data$Unseen)
                )
        )
}

#' Scale metrics
#'
#' Calculate baseline-scaled metric values
#'
#' @param metrics data.frame containing metrics
#' @param baseline_ranges date.frame containing baseline ranges
#'
#' @return data.frame containing summarised metrics
scale_metrics <- function(metrics, baseline_ranges) {
    metrics |>
        dplyr::left_join(
            baseline_ranges,
            by = c("Dataset", "Metric", "Type")
        ) |>
        dplyr::mutate(
            ScaledValue = (.data$Value - .data$Lower) /
                (.data$Upper - .data$Lower)
        )
}

#' Average method
#'
#' Compute the average metric scores for a given method
#'
#' @param metrics
#' @param method_pattern
#'
#' @details
#' Method names are matched to `method_pattern` using `stringr::str_detect()`
#' and any that match are replaced with `method_pattern`. The `data.frame` is
#' then group by `Dataset`, `Method`, `Integration`, `Type`, and `Metric` and
#' the mean value is calculated.
#'
#' @return data.frame containing the metric scores averaged for matching methods
average_method <- function(metrics, method_pattern) {
    metrics|>
        dplyr::mutate(
            Method = if_else(
                str_detect(.data$Method, method_pattern),
                method_pattern,
                .data$Method
            )
        ) |>
        dplyr::group_by(
            .data$Dataset,
            .data$Method,
            .data$Integration,
            .data$Type,
            .data$Metric
        ) |>
        dplyr::summarise(
            Value = mean(.data$Value),
            .groups = "drop"
        )
}
