#' BRCore ggplot2 Theme
#'
#' A consistent \pkg{ggplot2} theme used across all BRCore plotting functions.
#' Built on top of \code{\link[ggplot2]{theme_classic}} with standardised font
#' sizes and element dimensions.
#'
#' @param base_size Numeric. Base font size in points (default \code{11}).
#'   All other sizes are derived from this value.
#' @param axis_text_angle Numeric. Rotation angle for x-axis tick labels
#'   (default \code{0}). Set to e.g. \code{45} or \code{90} when labels overlap.
#' @param legend_position Character or numeric vector passed to
#'   \code{legend.position} (default \code{"right"}).
#' @param dashed_grid Logical. When \code{TRUE}, adds a dashed horizontal major
#'   grid line (used by rarefaction metric plots). Default \code{FALSE}.
#' @param extra_themes A named list of additional \pkg{ggplot2}
#'   \code{\link[ggplot2]{theme}} objects to compose onto the base theme.
#'   Useful for one-off per-plot overrides without breaking the shared theme.
#'   Example: \code{list(legend = theme(legend.title = element_blank()))}.
#'
#' @details
#' Size conventions derived from \code{base_size}:
#' \itemize{
#'   \item \strong{Plot title} – \code{base_size + 1}, bold, centred.
#'   \item \strong{Plot subtitle} – \code{base_size - 2}, centred.
#'   \item \strong{Axis titles} – \code{base_size}, bold.
#'   \item \strong{Axis text} – \code{base_size - 1}.
#'   \item \strong{Legend title} – \code{base_size - 1}, bold.
#'   \item \strong{Legend text} – \code{base_size - 2}.
#'   \item \strong{Strip text (facets)} – \code{base_size}, bold.
#'   \item \strong{Legend keys} – 0.4 cm × 0.4 cm.
#' }
#'
#' @return A \pkg{ggplot2} \code{\link[ggplot2]{theme}} object that can be added
#'   to any \code{ggplot} with \code{+}.
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   .brcore_theme()
#'
#' # With dashed grid and an extra per-plot override
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   .brcore_theme(
#'     dashed_grid  = TRUE,
#'     extra_themes = list(no_legend = theme(legend.position = "none"))
#'   )
#' }
#'
#' @importFrom ggplot2 theme_classic theme element_text element_line
#'   element_blank element_rect
#' @importFrom grid unit
#'
#' @noRd
#' @keywords internal
.brcore_theme <- function(
  base_size = 11,
  axis_text_angle = 0,
  legend_position = "right",
  dashed_grid = FALSE,
  extra_themes = list()
) {
  base <- theme_classic(base_size = base_size) +
    theme(
      # Titles ----
      plot.title = element_text(
        hjust = 0.5,
        size = base_size + 1,
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = base_size - 2,
        colour = "black"
      ),

      # Axes ----
      axis.title = element_text(
        size = base_size,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = base_size - 1,
        colour = "black"
      ),
      axis.text.x = element_text(
        angle = axis_text_angle,
        hjust = if (axis_text_angle == 0) 0.5 else 1,
        vjust = if (axis_text_angle == 0) 0.5 else 1,
        size = base_size - 1
      ),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),

      # Legend ----
      legend.position = legend_position,
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.key.width = grid::unit(0.4, "cm"),
      legend.title = element_text(
        size = base_size - 1,
        face = "bold",
        colour = "black"
      ),
      legend.text = element_text(
        size = base_size - 2,
        colour = "black"
      ),

      # Facet strips ----
      strip.background = element_blank(),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        colour = "black"
      ),

      # Panel ----
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 1)
    )

  # Optional dashed horizontal grid (e.g. rarefaction metric plots)
  if (dashed_grid) {
    base <- base +
      theme(
        panel.grid.major.y = element_line(
          color = "grey",
          linetype = "dashed",
          linewidth = 0.5
        )
      )
  }

  # Compose any extra per-plot theme overrides
  for (t in extra_themes) {
    base <- base + t
  }

  base
}
