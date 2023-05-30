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
                                     "IntegrationBatch",
                                     "IntegrationBio",
                                     "Mapping",
                                     "Classification",
                                     "Unseen"
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
        theme_report()
}

#' Plot methods dotplot
#'
#' Plot a dotplot showing the scores for different methods for a metric
#'
#' @param metrics Data frame containing metric scores
#' @param dataset Name of the dataset to plot scores for
#' @param metric Name of the metric to plot scores for
#'
#' @return ggplot2 object
plot_methods_dotplot <- function(metrics, dataset, metric) {

    metrics <- metrics |>
        dplyr::filter(
            .data$Dataset == dataset,
            .data$Metric == metric
        ) |>
        dplyr::mutate(
            IntegrationType = stringr::str_remove(
                .data$Integration,
                "-[0-9]+"
            )
        )

    ggplot2::ggplot(
        metrics,
        ggplot2::aes(
            x      = .data$Value,
            y      = .data$Method,
            colour = .data$IntegrationType
        )
    ) +
        ggplot2::geom_jitter(
            size   = 2,
            alpha  = 0.3,
            height = 0.2,
            width  = 0
        ) +
        ggplot2::stat_summary(
            fun.min = {\(x) {mean(x) - sd(x)}},
            fun     = mean,
            fun.max = {\(x) {mean(x) + sd(x)}},
            size    = 1,
            alpha   = 1,
            geom    = "linerange"
        ) +
        ggplot2::stat_summary(
            fun     = mean,
            size    = 3,
            alpha   = 1,
            geom    = "point"
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::facet_wrap(~ IntegrationType) +
        ggplot2::labs(title = metric) +
        theme_report(base_size = 12) +
        theme(
            legend.position = "none",
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()
        )
}

#' Plot methods variance
#'
#' Plot a bar chart showing the variance in scores for different methods for a
#' metric.
#'
#' @param metrics Data frame containing metric scores
#' @param dataset Name of the dataset to plot scores for
#' @param metric Name of the metric to plot scores for
#'
#' @details
#' This function assumes there are multiple values for each method (probably
#' from multiple integration runs). A bar for the between method variance is
#' also shown along with results from a one-way ANOVA.
#'
#' @return ggplot2 object
plot_methods_variance <- function(metrics, dataset, metric) {

    scores <- metrics |>
        dplyr::filter(
            .data$Dataset == dataset,
            .data$Metric == metric
        ) |>
        dplyr::mutate(
            IntegrationType = stringr::str_remove(
                .data$Integration,
                "-[0-9]+"
            )
        )

    summaries <- scores |>
        dplyr::group_by(Method, IntegrationType) |>
        dplyr::summarise(
            Mean     = mean(Value),
            Variance = var(Value),
            .groups = "drop"
        )

    grand_summaries <- summaries |>
        dplyr::group_by(IntegrationType) |>
        dplyr::summarise(
            Variance = var(Mean),
            Mean     = mean(Mean)
        ) |>
        dplyr::mutate(Method = "Between methods")

    summaries <- dplyr::bind_rows(summaries, grand_summaries) |>
        dplyr::mutate(
            Type = dplyr::if_else(
                Method == "Between methods",
                "Overall",
                "Single method"
            )
        ) |>
        dplyr::mutate(
            Method = forcats::fct_relevel(
                Method,
                "Between methods",
                after = Inf
            )
        )

    anovas <- scores |>
        dplyr::group_by(IntegrationType) |>
        dplyr::group_split() |>
        purrr::map_dfr(function(.data) {
            anova_res <- oneway.test(
                Value ~ Method,
                data      = .data,
                var.equal = TRUE
            )

            data.frame(
                IntegrationType = unique(.data$IntegrationType),
                Statistic       = anova_res$statistic,
                PValue          = anova_res$p.value
            )
        }) |>
        dplyr::mutate(
            Label = paste0(
                "One-way ANOVA\n",
                "F = ", signif(.data$Statistic, 5),
                ", p-value = ", signif(.data$PValue, 2)
            )
        )

    ggplot2::ggplot(
        summaries,
        ggplot2::aes(
            x      = .data$Variance,
            y      = .data$Method
        )
    ) +
        ggplot2::geom_col(
            ggplot2::aes(fill   = .data$Type)
        ) +
        ggplot2::geom_text(
            data = anovas,
            ggplot2::aes(label = .data$Label),
            x     = Inf, y     = -Inf,
            hjust = 1.1, vjust = -0.5
        ) +
        ggplot2::facet_wrap(~ IntegrationType) +
        ggplot2::labs(
            title = metric,
            x     = "Variance"
        ) +
        theme_report(base_size = 12) +
        theme(
            legend.position = "none",
            axis.title.y = ggplot2::element_blank()
        )
}

#' Theme report
#'
#' ggplot2 theme for report plots
#'
#' @param base_size Base font size (pts)
#' @param base_family Base font family
#' @param base_line_size Base size for line elements
#' @param base_rect_size Base size for rect elements
#' @param plot_title_size Size for plot title text
#' @param border Whether to add a border around the plot area
#' @param border_col Colour for border
#' @param grid Whether to add a grid to the plot, either `"both"` (all grids),
#' `"none"` (no grids), `"x"` (x-axis only) or `"y"` (y-axis only)
#' @param grid_col Colour for grid lines
#'
#' @details
#' Inspired by `hbrthemes::theme_ipsum()`
#'
#' @return theme object
theme_report <- function(base_size = 8, base_family = "",
                         base_line_size = base_size / 22,
                         base_rect_size = base_size / 22,
                         plot_title_size = base_size * 1.6,
                         border = TRUE, border_col = "grey30",
                         grid = c("both", "none", "x", "y"),
                         grid_col = "grey90") {

    grid <- match.arg(grid)

    theme <- ggplot2::theme_minimal(
        base_size, base_family, base_line_size, base_rect_size
    )

    if (border) {
        theme <- theme + ggplot2::theme(
            panel.border = ggplot2::element_rect(
                colour = border_col,
                fill   = "NA"
            )
        )
    }

    if (grid == "none") {
        theme <- theme + ggplot2::theme(panel.grid = ggplot2::element_blank())
    } else {

        theme <- theme + ggplot2::theme(
            panel.grid = ggplot2::element_line(
                colour = grid_col,
                size   = 0.3
            ),
            panel.grid.major = ggplot2::element_line(
                colour = grid_col,
                size   = 0.3
            ),
            panel.grid.minor = ggplot2::element_line(
                colour = grid_col,
                size   = ggplot2::rel(0.5)
            )
        )

        if (grid == "x") {
            theme <- theme + ggplot2::theme(
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank()
            )
        }

        if (grid == "y") {
            theme <- theme + ggplot2::theme(
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank()
            )
        }
    }

    theme + ggplot2::theme(
        plot.title = ggplot2::element_text(
            hjust  = 0,
            size   = plot_title_size,
            margin = ggplot2::margin(b = plot_title_size * 0.3),
            family = base_family,
            face   = "bold"
        ),
        axis.title.x = ggplot2::element_text(
            hjust  = 0,
            size   = base_size,
            family = base_family,
            face   = "bold"
        ),
        axis.title.y = ggplot2::element_text(
            hjust  = 0,
            size   = base_size,
            family = base_family,
            face   = "bold"
        ),
        legend.title = ggplot2::element_text(
            size   = base_size,
            family = base_family,
            face   = "bold",
            margin = ggplot2::margin(r = 10)
        ),
        strip.background = ggplot2::element_rect(fill = "black"),
        strip.text = element_text(
            size   = base_size,
            family = base_family,
            face   = "bold",
            colour = "white"
        )
    )
}
