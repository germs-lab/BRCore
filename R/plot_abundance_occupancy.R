#' Plot Abundance-Occupancy Curve and Display the Core Taxa
#'
#' Creates a scatter plot showing the relationship between mean relative abundance
#' and occupancy (occurrence frequency) of taxa, with core taxa highlighted.
#'
#' @param core_result A list object returned by \code{\link{identify_core}},
#'   containing at minimum:
#'   \itemize{
#'     \item \code{occupancy_abundance}: A data frame with columns \code{otu},
#'       \code{otu_rel} (mean relative abundance), and \code{otu_occ} (occupancy).
#'     \item \code{elbow_core}: Character vector of OTU IDs identified as core
#'       using the "elbow" method.
#'     \item \code{increase_core}: Character vector of OTU IDs identified as core
#'       using the "increase" method.
#'   }
#' @param core_set Character string specifying which core set to highlight.
#'   Must be either "elbow" or "increase" (Default elbow).
#'
#' @details The function creates a scatter plot where each point represents a
#'   taxon (i.e. ASV or OTU). Core taxa (as defined by the selected method) are
#'   shown in red, while non-core taxa are shown in grey. The plot uses a log10
#'   scale for abundance to better visualize the full range of abundances typically
#'    found in microbiome data.
#'
#' @return A ggplot object showing the abundance-occupancy plot with core taxa
#'   highlighted in red and non-core taxa in grey. The x-axis shows log10-transformed
#'   mean abundance and the y-axis shows occupancy (0-1).
#'
#' @seealso \code{\link{plot_core_distribution}} and \code{\link{identify_core}}
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' # Generate an object from `identify_core()` and then plot
#'
#' data("switchgrass", package = "BRCore")
#'
#' switchgrass_core <- identify_core(
#'   physeq_obj = switchgrass,
#'   priority_var = "sampling_date",
#'   increase_value = 0.02,
#'   abundance_weight = 0,
#'   seed = 1234
#' )
#'
#' plot_abundance_occupancy(core_result = switchgrass_core,
#'                          core_set = "increase")
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual theme_classic
#'   theme element_text element_rect element_blank labs alpha after_scale guides guide_legend
#' @importFrom grid unit
#'
#' @export
plot_abundance_occupancy <- function(core_result, core_set = "elbow") {
    # Validate core_set parameter
    if (!core_set %in% c("elbow", "increase")) {
        stop("core_set must be either 'elbow' or 'increase'")
    }

    # Select appropriate core OTU set
    if (core_set == "elbow") {
        core_otus <- core_result$elbow_core
    } else if (core_set == "increase") {
        core_otus <- core_result$increase_core
    }

    # Prepare data frame with membership classification
    whole_df <- core_result$abundance_occupancy |>
        mutate(
            membership = ifelse(otu %in% core_otus, "Core", "Not core")
        ) |>
        mutate(membership = as.factor(membership))

    # Create plot
    p <- whole_df |>
        ggplot(aes(
            x = log10(otu_rel),
            y = otu_occ,
            fill = membership
        )) +
        geom_point(
            shape = 21,
            size = 1.5,
            alpha = 0.9,
            aes(color = after_scale(fill))
        ) +
        scale_fill_manual(
            values = c("Core" = "#CC2D35", "Not core" = "grey")
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                size = 12,
                face = "bold"
            ),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            legend.position = c(0.98, 0.02),
            legend.justification = c("right", "bottom"),
            legend.background = element_rect(
                fill = alpha("white", 0.7),
                color = NA
            ),
            legend.key.height = unit(0.4, "cm"),
            legend.key.width = unit(0.4, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(size = 9)
        ) +
        labs(
            title = "Abundance-Occupancy curve",
            x = "Log10(mean abundance)",
            y = "Occupancy"
        )

    # change legend shape to avoid confusion with the scatter points.
    p <- p +
        guides(
            fill = guide_legend(
                override.aes = list(
                    shape = 22,
                    size = 5,
                    colour = NA,
                    stroke = 0
                )
            )
        )

    return(p)
}
