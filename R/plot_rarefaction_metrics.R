#' Plot pre-rarefaction diagnostics
#'
#' This function creates a 6-panel diagnostic plot showing sequencing depth,
#' Good's coverage, and outlier behavior based on a data frame or a
#' `phyloseq` object with rare stats already added via `add_rare_stats()`.
#'
#' @param data Either a `phyloseq` object or a `data.frame` that includes columns
#' `read_num`, `goods_cov`, and `outlier`.
#'
#' @return A `ggarrange` object with six plots.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point geom_jitter
#' @importFrom ggplot2 geom_boxplot geom_bar theme_bw theme labs
#' @importFrom ggplot2 scale_y_log10 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_classic element_line coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr ggarrange
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr filter arrange
#' @importFrom scales comma
#' @importFrom cli cli_alert_info cli_progress_step cli_alert_success cli_abort
#' @importFrom phyloseq sample_data

#'
#' @export
#'
plot_rarefaction_metrics <- function(data) {
    # Extract sample data depending on class
    if (inherits(data, "phyloseq")) {
        sample_df <- data.frame(sample_data(data))
        sample_df <- rownames_to_column(sample_df, "sample_id")
    } else if (inherits(data, "data.frame")) {
        sample_df <- data
    } else {
        stop("Input must be a phyloseq object or a data.frame")
    }

    # Force to tibble/data.frame
    sample_df <- as_tibble(sample_df)

    cli::cli_alert_info("Processing {nrow(sample_df)} sample{?s}")

    # Check required columns
    required_cols <- c("read_num", "goods_cov", "outlier")
    if (!all(required_cols %in% colnames(sample_df))) {
        cli::cli_abort(
            "Input data must contain columns: {.field read_num}, {.field goods_cov}, and {.field outlier}"
        )
    }

    # Calculate first quartile of read_num
    readno_q1 <- quantile(sample_df$read_num, probs = 0.25, na.rm = TRUE)

    cli::cli_progress_step("Generating rarefaction diagnostic plots")

    # Generate plot grid
    plot_output <- ggarrange(
        # Plot a - Histogram with improvements
        ggplot(sample_df, aes(x = read_num)) +
            geom_histogram(
                binwidth = 5000,
                fill = "firebrick",
                color = "white"
            ) +
            scale_x_continuous(labels = scales::comma) +
            scale_y_continuous(labels = scales::comma) +
            theme_classic() +
            theme(
                panel.grid.major.y = element_line(
                    color = "grey",
                    linetype = "dashed",
                    linewidth = 0.5
                )
            ) +
            labs(
                title = "Histogram",
                x = "Sequence reads",
                y = "Sample counts"
            ),

        # Plot b - Lower 25% histogram
        ggplot(sample_df, aes(x = read_num)) +
            geom_histogram(
                binwidth = 1000,
                fill = "firebrick",
                color = "white"
            ) +
            coord_cartesian(xlim = c(0, readno_q1)) +
            scale_x_continuous(labels = scales::comma) +
            scale_y_continuous(labels = scales::comma) +
            theme_classic() +
            theme(
                panel.grid.major.y = element_line(
                    color = "grey",
                    linetype = "dashed",
                    linewidth = 0.5
                )
            ) +
            labs(
                title = "Histogram (Lower 25%)",
                x = "Sequence reads",
                y = "Sample counts"
            ),

        # Plot c - Good's Coverage
        ggplot(sample_df, aes(x = read_num, y = goods_cov)) +
            geom_point(shape = 19, color = "firebrick", size = 1) +
            scale_x_continuous(labels = scales::comma) +
            theme_classic() +
            theme(
                panel.grid.major.y = element_line(
                    color = "grey",
                    linetype = "dashed",
                    linewidth = 0.5
                )
            ) +
            labs(
                title = "Good's Coverage",
                x = "Sequence reads",
                y = "Good's coverage %"
            ),

        # Plot d - Log10 jitter
        ggplot(sample_df, aes(x = 1, y = read_num)) +
            geom_jitter(
                shape = 19,
                color = "firebrick",
                width = 0.2,
                size = 1
            ) +
            theme_classic() +
            theme(
                panel.grid.major.y = element_line(
                    color = "grey",
                    linetype = "dashed",
                    linewidth = 0.5
                )
            ) +
            scale_y_log10(labels = scales::comma) +
            labs(title = "Log10 jitter", x = "Data set", y = "Sequence reads"),

        # Plot e - Log10 boxplot
        ggplot(sample_df, aes(x = 1, y = read_num)) +
            geom_boxplot(color = "firebrick", fill = "firebrick", alpha = 0.3) +
            geom_text_repel(
                data = filter(sample_df, !is.na(outlier)),
                mapping = aes(x = 1, y = read_num, label = outlier),
                max.overlaps = 15,
                size = 3
            ) +
            theme_classic() +
            theme(
                panel.grid.major.y = element_line(
                    color = "grey",
                    linetype = "dashed",
                    linewidth = 0.5
                )
            ) +
            scale_y_log10(labels = scales::comma) +
            labs(title = "Log10 boxplot", x = "Data set", y = "Sequence reads"),

        # Plot f - Ranked samples (fixed)
        {
            temp_df <- arrange(sample_df, read_num)
            ggplot(temp_df, aes(x = 1:nrow(temp_df), y = read_num)) +
                geom_bar(stat = "identity", fill = "firebrick", color = NA) +
                scale_y_continuous(labels = scales::comma) +
                theme_classic() +
                theme(
                    panel.grid.major.y = element_line(
                        color = "grey",
                        linetype = "dashed",
                        linewidth = 0.5
                    )
                ) +
                labs(
                    title = "Ranked samples",
                    x = "Samples",
                    y = "Sequence reads"
                )
        },

        ncol = 3,
        nrow = 2,
        align = "hv",
        labels = c("a", "b", "c", "d", "e", "f")
    )

    # Display success message
    cli::cli_alert_success(
        "Rarefaction diagnostic plots generated successfully"
    )

    return(plot_output)
}
