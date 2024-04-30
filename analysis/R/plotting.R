#' Plot metrics counts
#'
#' Plot a heatmap of the counts of metrics by integration, dataset, and method
#'
#' @param metrics data.frame containing metrics
#'
#' @return A `patchwork` object
plot_metric_counts <- function(metrics) {

    n_integrations <- length(unique(metrics$Integration))
    n_datasets <- length(unique(metrics$Dataset))
    n_methods <- length(unique(metrics$Method))
    n_total <- n_integrations + n_datasets + n_methods

    counts_integration <- metrics |>
        dplyr::group_by(.data$Integration, .data$Metric) |>
        dplyr::count(name = "Count")

    counts_dataset <- metrics |>
        dplyr::group_by(.data$Dataset, .data$Metric) |>
        dplyr::count(name = "Count")

    counts_method <- metrics |>
        dplyr::group_by(.data$Method, .data$Metric) |>
        dplyr::count(name = "Count")

    counts_integration_plot <- ggplot2::ggplot(
        counts_integration,
        ggplot2::aes(
            x = .data$Integration, y = forcats::fct_rev(.data$Metric),
            fill = .data$Count
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By integration"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom"
        )

    counts_dataset_plot <- ggplot2::ggplot(
        counts_dataset,
        ggplot2::aes(
            x = .data$Dataset, y = forcats::fct_rev(.data$Metric),
            fill = .data$Count
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By dataset"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )

    counts_method_plot <- ggplot2::ggplot(
        counts_method,
        ggplot2::aes(
            x = .data$Method, y = forcats::fct_rev(.data$Metric),
            fill = .data$Count
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By method"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )

    patchwork::wrap_plots(
        counts_integration_plot, counts_dataset_plot, counts_method_plot,
        widths = c(n_integrations, n_datasets, n_methods) / n_total
    )
}

#' Plot metrics missing
#'
#' Plot a heatmap of the count of missing values of metrics by integration,
#' dataset, and method
#'
#' @param metrics data.frame containing metrics
#'
#' @return A `patchwork` object
plot_metric_missing <- function(metrics) {

    n_integrations <- length(unique(metrics$Integration))
    n_datasets <- length(unique(metrics$Dataset))
    n_methods <- length(unique(metrics$Method))
    n_total <- n_integrations + n_datasets + n_methods

    missing_integration <- metrics |>
        dplyr::group_by(.data$Integration, .data$Metric) |>
        dplyr::summarise(
            Missing = sum(is.na(.data$Value)),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            Missing = dplyr::if_else(.data$Missing == 0, NA, .data$Missing)
        )

    missing_dataset <- metrics |>
        dplyr::group_by(.data$Dataset, .data$Metric) |>
        dplyr::summarise(
            Missing = sum(is.na(.data$Value)),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            Missing = dplyr::if_else(.data$Missing == 0, NA, .data$Missing)
        )

    missing_method <- metrics |>
        dplyr::group_by(.data$Method, .data$Metric) |>
        dplyr::summarise(
            Missing = sum(is.na(.data$Value)),
            .groups = "drop"
        ) |>
        mutate(
            Missing = dplyr::if_else(.data$Missing == 0, NA, .data$Missing)
        )

    missing_integration_plot <- ggplot2::ggplot(
        missing_integration,
        ggplot2::aes(
            x = .data$Integration, y = forcats::fct_rev(.data$Metric),
            fill = .data$Missing
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By integration"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom"
        )

    missing_dataset_plot <- ggplot2::ggplot(
        missing_dataset,
        ggplot2::aes(
            x = .data$Dataset, y = forcats::fct_rev(.data$Metric),
            fill = .data$Missing
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By dataset"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )

    missing_method_plot <- ggplot2::ggplot(
        missing_method,
        ggplot2::aes(
            x = .data$Method, y = forcats::fct_rev(.data$Metric),
            fill = .data$Missing
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(
            title = "By method"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )

    patchwork::wrap_plots(
        missing_integration_plot, missing_dataset_plot, missing_method_plot,
        widths = c(n_integrations, n_datasets, n_methods) / n_total
    )
}

#' Theme features
#'
#' ggplot2 theme for feature selection analysis
#'
#' @param ... Additional arguments passed `ggplot2::theme_minimal()`
#'
#' @return **{ggplot2}** `theme` object
theme_features <- function(...) {
    ggplot2::theme_minimal(...) +
        ggplot2::theme(
            panel.border = ggplot2::element_rect(fill = NA),
            strip.text = ggplot2::element_text(colour = "white"),
            strip.background = ggplot2::element_rect(fill = "black")
        )
}

#' Theme features pub
#'
#' ggplot2 theme for feature selection publication figures with set text sizes
#'
#' @param ... Additional arguments passed `theme_features()`
#'
#' @return **{ggplot2}** `theme` object
theme_features_pub <- function(...) {
    theme_features(base_size = 8, ...) +
        theme(
            plot.title = element_text(size = 8),
            strip.text = element_text(size = 6),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 5),
            legend.key.size = unit(0.4, "cm")
        )
}

#' Save figure files
#'
#' Save output figure files
#'
#' @param figure `ggplot` object containing figure to save
#' @param base_path Base path to save files to. File extensions to be appended
#' using `fs::path_ext_set()`
#' @param width Width of the output figure. Cannot exceed A4 portrait width.
#' @param height Height of the output figure.
#'
#' @return `base_path` invisibly
save_figure_files <- function(figure, base_path, width = 8, height = 6) {
    if (width > 8.3) {
        stop("Width cannot be greater than A4 portrait width (8.3 inches)")
    }

    png_path <- fs::path_ext_set(base_path, "png")
    ggsave(png_path, figure, width = width, height = height)

    pdf_path <- fs::path_ext_set(base_path, "pdf")
    ggsave(pdf_path, figure, width = width, height = height)

    invisible(base_path)
}
