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
#' @param dataset_name Optional character string. When provided, it is
#'   prepended to the plot title (e.g. \code{"Switchgrass"}).
#'   Default \code{NULL} (no prefix).
#'
#' @return A list containing: 1) `df_for_plot`, a data frame used for plotting, and 2) `plot_identified_core`, a ggplot object visualizing the Bray-Curtis increase with annotated cutoff points.
#'
#' @details
#' The function converts `rank` to integers and zooms the x-axis to the first
#' `1.2 * lastCall` ranks. Label positions are computed dynamically from the
#' observed `proportionBC` range to avoid overlap.
#'
#' @seealso \link[=identify_core]{identify_core()}
#'
#' @examples
#' library(BRCore)
#' data("switchgrass_core", package = "BRCore")
#'
#' p <- plot_identified_core(
#'   bray_curtis_ranked = switchgrass_core$bray_curtis_ranked,
#'   elbow = switchgrass_core$elbow,
#'   lastCall = switchgrass_core$bc_increase,
#'   increase_value = switchgrass_core$increase_value
#' )
#' print(p)
#'
#' @seealso \code{\link{identify_core}}
#'
#' @importFrom grid unit
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot aes geom_point geom_vline annotate labs
#' @importFrom ggplot2 theme_classic theme element_text scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous expansion
#'
#' @export
plot_identified_core <- function(
  bray_curtis_ranked,
  elbow,
  lastCall,
  increase_value = 0.02,
  dataset_name = NULL
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
    labs(
      title = if (!is.null(dataset_name)) {
        bquote(bold(
          "Core Microbiome Identification for" ~ italic(
            .(dataset_name)
          ) ~
            "dataset"
        ))
      } else {
        expression(bold("Core Microbiome Identification"))
      },
      subtitle = paste0(
        "Identification Methods: Elbow (",
        elbow,
        " OTUs), Last ",
        percent_label,
        " (",
        lastCall,
        " OTUs)"
      ),
      x = "Ranked OTUs",
      y = "% Bray-Curtis similarity"
    ) +
    .brcore_theme() +
    scale_x_continuous(breaks = seq(0, max(new_data$rank_num), by = 5))

  return(list(
    df_for_plot = new_data,
    plot_identified_core = p
  ))
}
