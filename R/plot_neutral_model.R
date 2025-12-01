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
#' 
#' The function does not recompute statistics; it only visualizes the supplied
#' predictions and metrics.
#' 
#' @examples
#' \dontrun{
#' # Assuming you have run a neutral model fit on the switchgrass data.frames
#' is.data.frame(switchgrass_core_fit$model_prediction)
#' is.data.frame(switchgrass_core_fit$goodness_of_fit)
#' 
#' p <- plot_neutral_model(switchgrass_core_fit)
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
    obs1 <- as.data.frame(fit_result$model_prediction)
    obs2 <- as.data.frame(fit_result$goodness_of_fit)
    
    # pre-split - now handling all 6 combinations
    pts_notcore_as    <- dplyr::filter(obs1, membership == "Not core", fit_class == "As predicted")
    pts_notcore_above <- dplyr::filter(obs1, membership == "Not core", fit_class == "Above prediction")
    pts_notcore_below <- dplyr::filter(obs1, membership == "Not core", fit_class == "Below prediction")
    pts_core_as       <- dplyr::filter(obs1, membership == "Core", fit_class == "As predicted")
    pts_core_above    <- dplyr::filter(obs1, membership == "Core", fit_class == "Above prediction")
    pts_core_below    <- dplyr::filter(obs1, membership == "Core", fit_class == "Below prediction")
    
    # safe combine - all "Not core" groups get the same label for legend grouping
    all_pts <- dplyr::bind_rows(
        dplyr::mutate(pts_notcore_as,    grp = "Not core"),
        dplyr::mutate(pts_notcore_above, grp = "Not core"),
        dplyr::mutate(pts_notcore_below, grp = "Not core"),
        dplyr::mutate(pts_core_as,       grp = "Core (as predicted)"),
        dplyr::mutate(pts_core_above,    grp = "Core (above prediction)"),
        dplyr::mutate(pts_core_below,    grp = "Core (below prediction)")
    )
    
    # if still empty, bail nicely
    if (nrow(all_pts) == 0) {
        cli::cli_abort("Nothing to plot: no rows in model_prediction after filtering (check p>0 and membership/fit_class).")
    }
    
    #print(all_pts)
    
    # Extract R2 and m (works for column- or rowname-shaped tables)
    R2 <- if ("Rsqr" %in% names(obs2)) obs2$Rsqr[1] else suppressWarnings(as.numeric(obs2["Rsqr", 1]))
    mV <- if ("m"    %in% names(obs2)) obs2$m[1]    else suppressWarnings(as.numeric(obs2["m",    1]))
    
    # Compute a spot just ABOVE the legend (bottom-right)
    xr <- range(log10(all_pts$otu_rel), na.rm = TRUE); dx <- diff(xr)
    yr <- range(all_pts$otu_occ,        na.rm = TRUE); dy <- diff(yr)
    box_x <- xr[1] + 0.95 * dx
    box_y <- yr[1] + 0.12 * dy   # raise/lower if you want different spacing
    
    p <- ggplot2::ggplot(all_pts, ggplot2::aes(x = log10(otu_rel), y = otu_occ, fill = grp)) +
        ggplot2::geom_point(shape = 21, size = 1.5, alpha = 0.9,
                            ggplot2::aes(color = after_scale(fill))) +
        ggplot2::scale_fill_manual(
            values = c(
                "Not core"                 = "grey",
                "Core (as predicted)"      = "#CC2D35",
                "Core (above prediction)"  = "#0072B2",
                "Core (below prediction)"  = "#7FB800"
            ),
            breaks = c("Not core", "Core (as predicted)",
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
            legend.title         = ggplot2::element_blank(),
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
            #label.size = 0  # This removes the border
        ) +
        ggplot2::labs(title = "Neutral model",
                      x = "Log10(mean abundance)",
                      y = "Occupancy")
    # change legend shape to avoid confusion with the scatter points.
    p <- p + guides(
        fill = guide_legend(
            override.aes = list(shape = 22, size = 4, colour = NA, stroke = 0)
        )
    )
    
    return(p)
    
}

