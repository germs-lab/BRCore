#' Plot Bray-Curtis increase over ranked OTU/ASVs
#'
#' @description
#' Visualize the cumulative normalized mean Bray-Curtis increase returned by
#' \link[=identify_core]{identify_core()}, over ranked OTU/ASVs and shows cutoff
#' points for elbow percent increase methods.
#'
#' @param bray_curtis_ranked A tibble as returned by `identify_core()$bray_curtis_ranked`.
#' @param elbow The number of OTU/ASVs identified by the elbow method (Integer).
#' @param lastCall The number of OTU/ASVs identified by the last percent
#' Bray-Curtis increase method (Integer).
#' @param increase_value The percent increase value in decimal (e.g. 0.02) used
#' for the Bray-Curtis increase method.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @details
#' The function converts `rank` to integers and zooms the x-axis to the first
#' `1.2 * lastCall` ranks. Label positions are computed dynamically from the
#' observed `proportionBC` range to avoid overlap.
#'
#' @seealso \link[=identify_core]{identify_core()}
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' # Example with the package switchgrass dataset
#' data("switchgrass", package = "BRCore")
#'
#' # Identify core taxa
#' res <- identify_core(
#'  physeq_obj = switchgrass,
#'  priority_var = "sampling_date",
#'  increase_value = 0.02,
#'  seed = 48821
#' )
#'
#' # Plot using the returned curve and cut indices; label from increase_value
#' plot_identified_core(
#'   bray_curtis_ranked = res$bray_curtis_ranked,
#'   elbow = res$elbow,
#'   lastCall = res$bc_increase,
#'   increase_value = res$increase_value
#'  )
#' }
#'
#' @seealso \code{\link{identify_core}}
#'
#' @importFrom grid unit
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot aes geom_point geom_vline annotate labs theme_classic theme element_text scale_x_continuous scale_y_continuous expansion
#' @export

plot_identified_core <- function(
    bray_curtis_ranked,
    elbow,
    lastCall,
    increase_value = 0.02
) {
    new_data <- bray_curtis_ranked |>
        mutate(
            rank_num = suppressWarnings(as.integer(as.character(rank)))
        ) |>
        filter(rank_num <= 1.2 * lastCall) |>
        mutate(rank_fac = factor(rank_num, levels = rank_num))

    # Calculate percentage label text
    percent_label <- paste0(round(increase_value * 100), "%")

    # Calculate dynamic label positions
    max_rank <- max(new_data$rank_num, na.rm = TRUE)
    max_prop <- max(new_data$proportionBC, na.rm = TRUE)
    min_prop <- min(new_data$proportionBC, na.rm = TRUE)
    prop_range <- max_prop - min_prop

    # Dynamic y positions for labels
    elbow_y <- min_prop + 0.15 * prop_range
    lastCall_y <- min_prop + 0.35 * prop_range

    # Adjust x position for labels if too close to edges
    elbow_x_offset <- ifelse(
        elbow < max_rank * 0.1,
        elbow + max_rank * 0.02,
        elbow - max_rank * 0.02
    )
    lastCall_x_offset <- ifelse(
        lastCall > max_rank * 0.9,
        lastCall - max_rank * 0.02,
        lastCall + max_rank * 0.02
    )

    # Base plot
    p <- ggplot(new_data, aes(x = rank_num, y = proportionBC)) +
        geom_point(size = 1.5, alpha = 0.7) +

        # Enhanced vertical lines with better visibility
        geom_vline(
            xintercept = elbow,
            linetype = "dashed",
            color = "#d73027",
            linewidth = 1.2,
            alpha = 0.8
        ) +
        geom_vline(
            xintercept = lastCall,
            linetype = "dashed",
            color = "#1a9850",
            linewidth = 1.2,
            alpha = 0.8
        ) +

        # # Background label rectangles and positioning
        annotate(
            "rect",
            xmin = 0,
            xmax = elbow,
            ymin = -Inf,
            ymax = Inf,
            fill = "white",
            alpha = 0.05
        ) +
        annotate(
            "rect",
            xmin = elbow,
            xmax = lastCall,
            ymin = -Inf,
            ymax = Inf,
            fill = "white",
            alpha = 0.05
        ) +

        annotate(
            "label",
            x = elbow_x_offset,
            y = elbow_y,
            label = paste0("Elbow (", elbow, ")"),
            color = "#d73027",
            fill = "white",
            size = 3.5,
            fontface = "bold",
            label.padding = unit(0.3, "lines")
        ) +
        annotate(
            "label",
            x = lastCall_x_offset,
            y = lastCall_y,
            label = paste0("Last ", percent_label, " (", lastCall, ")"),
            color = "#1a9850",
            fill = "white",
            size = 3.5,
            fontface = "bold",
            label.padding = unit(0.3, "lines")
        ) +

        # Titles and labels
        labs(
            title = paste0(
                "Core Microbiome Identification: Elbow and Last ",
                percent_label,
                " Methods"
            ),
            subtitle = paste0(
                "Core taxa identified: Elbow method (",
                elbow,
                " OTUs), Last ",
                percent_label,
                " method (",
                lastCall,
                " OTUs)"
            ),
            x = "Ranked OTUs",
            y = "% Bray-Curtis similarity"
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9)
        ) +
        # Scaling
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

    # To add some more fancy graphics later...
    # axis.title = element_text(size = 11),
    # axis.text = element_text(size = 9),
    # panel.grid.minor = element_blank(),
    # panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    # panel.grid.major.y = element_line(color = "gray90", size = 0.3))

    return(p)
}
