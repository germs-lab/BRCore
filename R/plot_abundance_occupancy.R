#' Plot Abundance-Occupancy Relationship for Core Taxa
#'
#' Creates a scatter plot showing the relationship between mean relative abundance
#' and occupancy (occurrence frequency) of taxa, with core taxa highlighted.
#'
#' @param core_result A list object returned by a core identification function,
#'   containing at minimum:
#'   \itemize{
#'     \item \code{occupancy_abundance}: A data frame with columns \code{otu},
#'       \code{otu_rel} (mean relative abundance), and \code{otu_occ} (occupancy).
#'     \item \code{elbow_core}: Character vector of OTU IDs identified as core
#'       using the elbow method.
#'     \item \code{increase_core}: Character vector of OTU IDs identified as core
#'       using the increase method.
#'   }
#' @param core_set Character string specifying which core set to highlight.
#'   Must be either "elbow" or "increase" (Default elbow).
#'
#' @return A ggplot object showing the abundance-occupancy plot with core taxa
#'   highlighted in red and non-core taxa in grey. The x-axis shows log10-transformed
#'   mean abundance and the y-axis shows occupancy (0-1).
#'
#' @details The function creates a scatter plot where each point represents a taxon.
#'   Core taxa (as defined by the selected method) are shown in red, while non-core
#'   taxa are shown in grey. The plot uses a log10 scale for abundance to better
#'   visualize the full range of abundances typically found in microbiome data.
#'
#' @examples
#' \dontrun{
#' # Generate an object from identify_core and then plot
#' 
#' data("switchgrass", package = "BRCore")
#' 
#' switchgrass_core <- identify_core(
#'   physeq_obj = switchgrass,
#'   priority_var = "sampling_date",  
#'   increase_value = 0.02,
#'   abundance_weight = 0,
#'   seed = 1234)
#' )
#'   
#' plot_abundance_occupancy(core_result =switchgrass_core,
#'                          core_set = "increase")
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual theme_classic
#'   theme element_text element_rect element_blank labs alpha after_scale
#' @importFrom grid unit
#' @export
plot_abundance_occupancy <- function(core_result, 
                                     core_set = "elbow") {
    
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
    whole_df <- core_result$occupancy_abundance %>% 
        dplyr::mutate(membership = base::ifelse(otu %in% core_otus, "Core", "Not core")) %>% 
        dplyr::mutate(membership = base::as.factor(membership))
    
    # Create plot
    whole_df %>% 
        ggplot2::ggplot(ggplot2::aes(x = base::log10(otu_rel), y = otu_occ, fill = membership)) +
        ggplot2::geom_point(shape = 21, size = 1.5, alpha = 0.9, 
                            ggplot2::aes(color = ggplot2::after_scale(fill))) +
        ggplot2::scale_fill_manual(values = c("Core" = "#CC2D35", "Not core" = "grey")) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
            legend.position = c(0.98, 0.02),
            legend.justification = c("right", "bottom"),
            legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.7), 
                                                      color = NA),
            legend.key.height = grid::unit(0.2, "cm"),
            legend.key.width = grid::unit(0.3, "cm"),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 8)
        ) +
        ggplot2::labs(
            title = "Abundance-Occupancy curve",
            x = "Log10(mean abundance)",
            y = "Occupancy"
        )
}