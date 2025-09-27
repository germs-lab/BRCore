#' Plot a fitted Neutral Model to Microbial Community Data
#'
#' This function plots ASV/OTUs log abundances into a fitted neutral model of 
#' microbial abundance-occupancy distribution. ASV/OTUs that are *Core*, *As predicted* (i.e. Neutral), 
#' *Below* and *Above* model predictions are drawn with distinct point colors, see details.
#' 
#' @param fit_result A list-like object returned by \link[=fit_neutral_model]{fit_neutral_model()}.
#' 
#' @return A \link[ggplot2:ggplot]{ggplot2} object.
#' 
#' @details
#' Points are split into four groups for display:
#' \itemize{
#'   \item \strong{Not core (as predicted)} – background points;
#'   \item \strong{Core (as predicted)} – core taxa whose occupancy matches the model;
#'   \item \strong{Core (above prediction)} – core taxa above the 95\% CI;
#'   \item \strong{Core (below prediction)} – core taxa below the 95\% CI.
#' }
#' 
#' @return
#' A \link[ggplot2:ggplot]{ggplot} object with:
#' \itemize{
#'   \item \emph{x}-axis: \code{log10(mean abundance)} (\code{log10(otu_rel)});
#'   \item \emph{y}-axis: \code{Occupancy} (\code{otu_occ});
#'   \item Points colored by membership/fit class as described above;
#'   \item Neutral-model curve (solid) with 95\% CI (grey, twodash);
#'   \item An inset white label reporting \emph{R}\eqn{^2} and \emph{m} taken
#'         directly from \code{fit_result$goodness_of_fit}.
#' }
#' The function does not recompute statistics; it only visualizes the supplied
#' predictions and metrics.
#' 
#' @examples
#' \dontrun{
#' # Assuming you have run a neutral model fit:
#' # fit_result$model_prediction and fit_result$goodness_of_fit are data.frames
#' p <- plot_neutral_model(fit_result)
#' print(p)
#' }
#' 
#' @seealso \code{\link{fit_neutral_model}}
#' 
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_classic theme
#'   element_text element_blank annotate 
#' @importFrom grid unit
#' @export
plot_neutral_model <- function(fit_result){
    
    # Export dataframes for plotting and annotations
    obs1 <- as.data.frame(fit_result$model_prediction) %>%
        dplyr::filter(!is.na(p))
    
    obs2 <- as.data.frame(fit_result$goodness_of_fit)
    
    # Extract R2 and m (works for column- or rowname-shaped tables)
    R2 <- if ("Rsqr" %in% names(obs2)) obs2$Rsqr[1] else suppressWarnings(as.numeric(obs2["Rsqr", 1]))
    mV <- if ("m"    %in% names(obs2)) obs2$m[1]    else suppressWarnings(as.numeric(obs2["m",    1]))
    
    # Pre-split points
    pts_bg         <- obs1[is.na(obs1$membership) | obs1$fit_class == "As predicted", , drop = FALSE]
    pts_core_as    <- obs1[obs1$membership == "core" & obs1$fit_class == "As predicted", , drop = FALSE]
    pts_core_above <- obs1[obs1$membership == "core" & obs1$fit_class == "Above prediction", , drop = FALSE]
    pts_core_below <- obs1[obs1$membership == "core" & obs1$fit_class == "Below prediction", , drop = FALSE]
    
    all_pts <- rbind(
        transform(pts_bg,         grp = "Not core (as predicted)"),
        transform(pts_core_as,    grp = "Core (as predicted)"),
        transform(pts_core_above, grp = "Core (above prediction)"),
        transform(pts_core_below, grp = "Core (below prediction)")
    )
    
    # Compute a spot just ABOVE the legend (bottom-right)
    xr <- range(log10(all_pts$otu_rel), na.rm = TRUE); dx <- diff(xr)
    yr <- range(all_pts$otu_occ,        na.rm = TRUE); dy <- diff(yr)
    box_x <- xr[1] + 0.95 * dx
    box_y <- yr[1] + 0.12 * dy   # raise/lower if you want different spacing
    
    ggplot2::ggplot(all_pts, ggplot2::aes(x = log10(otu_rel), y = otu_occ, fill = grp)) +
        ggplot2::geom_point(shape = 21, size = 1.5, alpha = 0.9,
                            ggplot2::aes(color = after_scale(fill))) +
        ggplot2::scale_fill_manual(
            values = c(
                "Not core (as predicted)"  = "grey",
                "Core (as predicted)"      = "#CC2D35",
                "Core (above prediction)"  = "#0072B2",
                "Core (below prediction)"  = "#7FB800"
            ),
            breaks = c("Not core (as predicted)", "Core (as predicted)",
                       "Core (above prediction)", "Core (below prediction)")
        ) +
        # Keep lines but remove their legends
        ggplot2::geom_line(data = obs1, ggplot2::aes(x = log10(p), y = freq.pred),
                           color = "blue", linewidth = 0.8, alpha = 0.25,
                           inherit.aes = FALSE, show.legend = FALSE) +
        ggplot2::geom_line(data = obs1, ggplot2::aes(x = log10(p), y = pred.upr),
                           color = "black", linewidth = 0.8, linetype = "twodash", alpha = 0.25,
                           inherit.aes = FALSE, show.legend = FALSE) +
        ggplot2::geom_line(data = obs1, ggplot2::aes(x = log10(p), y = pred.lwr),
                           color = "black", linewidth = 0.8, linetype = "twodash", alpha = 0.25,
                           inherit.aes = FALSE, show.legend = FALSE) +
        # Make sure only the Membership legend remains
        ggplot2::guides(color = "none", linetype = "none") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
            legend.position      = c(0.98, 0.02),
            legend.justification = c("right","bottom"),
            legend.background    = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.7), color = NA),
            legend.key.height    = grid::unit(0.2, "cm"),
            legend.key.width     = grid::unit(0.3, "cm"),
            legend.title = ggplot2::element_blank(),
            legend.text          = ggplot2::element_text(size = 8)
        ) +
        # White box with two lines: R^2 on top, m below (above the legend)
        ggplot2::annotate(
            "label", x = box_x, y = box_y,
            label = paste0(
                "atop(italic(R)^2==", sprintf("%.2f", R2), 
                ",italic(m)==",       sprintf("%.3f", mV), ")"
            ),
            parse = TRUE, hjust = 1, vjust = 0, size = 3,
            fill = "white", alpha = 0.9, lineheight = 1.05
        ) +
        ggplot2::labs(title = "Neutral model",
                      x = "Log10(mean abundance)",
                      y = "Occupancy")
}
