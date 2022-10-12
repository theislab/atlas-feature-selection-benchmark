#' Plot metrics dotplot
#'
#' Plot a dotplot of all metrics for a dataset
#'
#' @param metrics data.frame containing the metrics to plot
#' @param dataset Name of the dataset to plot
#' @param metric_type Type of metrics to plot
#' @param value Name of the column containing the values to plot
#'
#' @return ggplot object
plot_metrics_dotplot <- function(metrics, dataset,
                                 metric_type = c(
                                     "All",
                                     "Classification",
                                     "Integration"
                                 ),
                                 value = "Value"
                                ) {

    metric_type <- match.arg(metric_type)

    metrics <- dplyr::filter(metrics, .data$Dataset == dataset)

    if (metric_type != "All") {
        metrics <- dplyr::filter(metrics, .data$Type == metric_type)
    }

    ggplot2::ggplot(
        metrics,
        ggplot2::aes(
            x      = .data$Metric,
            y      = .data[[value]],
            colour = .data$Method,
            shape  = .data$Integration
        )
    ) +
        ggplot2::geom_jitter(
            size   = 5,
            width  = 0.3,
            height = 0,
        ) +
        # ggplot2::scale_shape_manual(
        #     values = c("circle filled", "triangle filled")
        # ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(shape = "circle filled")
            )
        ) +
        ggplot2::theme_minimal()
}
