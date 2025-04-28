#' Visualize Ordination Results with Flexible Components
#'
#' Plot ordinations reulst from dbRDA, capscale, PCoA and NMDS analyses.
#'
#' @param ordination_result Result object from vegan (dbRDA, capscale, rda, etc.)
#' @param metadata_df Dataframe containing sample metadata (must match ordination samples)
#' @param sample_id_col Column name in metadata containing sample IDs (default: "sample_id")
#' @param color_var Optional variable for point coloring
#' @param shape_var Optional variable for point shapes
#' @param arrow_scaling Scaling factor for constraint arrows (default: 2)
#' @param ellipse Logical - whether to draw ellipses (default: TRUE)
#' @param biplot Logical - whether to draw biplot arrows (default: TRUE)
#' @param biplot_labels Logical - whether to label biplot arrows (default: TRUE)
#' @param biplot_n Number of top constraints to show (default: 5)
#'
#' @section Acknowledgments:
#' This function was optimized with AI assistance from DeepSeek R1
#' for error handling, and performance improvements through iterative development and debugging (March 2025).
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' # For dbRDA:
#' brc_flex_ordi(model, metadata, sample_id_col = "sample_id")
#'
#' # For PCoA/NMDS:
#' brc_flex_ordi(ord_result, metadata, color_var = "Treatment")
#'
brc_flex_ordi <- function(
  ordination_result,
  metadata_df,
  sample_id_col = "sample_id",
  color_var = NULL,
  shape_var = NULL,
  arrow_scaling = 2,
  ellipse = TRUE,
  biplot = TRUE,
  biplot_labels = TRUE,
  biplot_n = 5
) {
  # Input validation
  if (!inherits(ordination_result, c("rda", "capscale", "metaMDS", "pcoa"))) {
    stop(
      "Input must be a vegan ordination object (rda, capscale, metaMDS, or pcoa)"
    )
  }

  # Extract scores
  site_scores <- as.data.frame(vegan::scores(
    ordination_result,
    display = "sites"
  ))
  site_scores[[sample_id_col]] <- rownames(site_scores)

  # Merge with metadata
  plot_data <- dplyr::left_join(
    site_scores,
    metadata_df,
    by = sample_id_col
  )

  # Determine axis labels
  if (inherits(ordination_result, c("rda", "capscale"))) {
    # dbRDA/RDA case
    summ <- summary(ordination_result)
    xlab <- paste0(
      "dbRDA1 (",
      round(summ$concont$importance[2, 1] * 100, 1),
      "%)"
    )
    ylab <- paste0(
      "dbRDA2 (",
      round(summ$concont$importance[2, 2] * 100, 1),
      "%)"
    )
  }
  if (inherits(ordination_result, "metaMDS")) {
    # NMDS case
    xlab <- "NMDS1"
    ylab <- "NMDS2"
  }
  if (inherits(ordination_result, "wcmdscale")) {
    # PCoA/default case
    xlab <- "PCo1"
    ylab <- "PCo2"
  }

  # Base plot
  p <- ggplot(
    plot_data,
    aes_string(
      x = colnames(site_scores)[1],
      y = colnames(site_scores)[2]
    )
  ) +
    geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey30", linetype = "dashed") +
    theme_classic(base_size = 12) +
    labs(x = xlab, y = ylab) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  # Add points with optional aesthetics
  point_aes <- aes()
  if (!is.null(color_var)) point_aes <- aes_string(color = color_var)
  if (!is.null(shape_var))
    point_aes <- modifyList(point_aes, aes_string(shape = shape_var))

  p <- p + geom_point(point_aes, size = 3, alpha = 0.7)

  # Add ellipses if requested
  if (ellipse && !is.null(color_var)) {
    p <- p +
      stat_ellipse(
        aes_string(color = color_var),
        geom = "path",
        linewidth = 1,
        type = "t",
        level = 0.95
      )
  }

  # Add biplot arrows if available and requested
  if (biplot && inherits(ordination_result, c("rda", "capscale"))) {
    arrow_scores <- as.data.frame(vegan::scores(
      ordination_result,
      display = "bp"
    ))

    # Select top constraints if requested
    if (!is.null(biplot_n) && nrow(arrow_scores) > biplot_n) {
      constraint_importance <- apply(
        arrow_scores[, 1:2],
        1,
        function(x) sqrt(sum(x^2))
      )
      arrow_scores <- arrow_scores[order(-constraint_importance), ][
        1:biplot_n,
      ]
    }

    arrow_scores <- arrow_scores * arrow_scaling

    p <- p +
      geom_segment(
        aes(x = 0, y = 0, xend = arrow_scores[, 1], yend = arrow_scores[, 2]),
        data = arrow_scores,
        arrow = arrow(length = unit(0.02, "npc")),
        color = "black",
        linewidth = 0.8
      )

    if (biplot_labels) {
      p <- p +
        geom_text(
          aes(
            x = arrow_scores[, 1] * 1.1,
            y = arrow_scores[, 2] * 1.1,
            label = rownames(arrow_scores)
          ),
          data = arrow_scores,
          size = 4,
          fontface = "bold"
        )
    }
  }

  return(p)
}
