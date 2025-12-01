#' Plot Core Taxa  Occupancy Across Metadata Groups
#'
#' Creates a bar plot showing core taxa (i.e. OTUs/ASVs) occupancy patterns across 
#' a grouping variable.
#' 
#' @param core_result A list object returned by \code{\link{identify_core}},
#'   containing at minimum:
#'   \itemize{
#'     \item \code{otu_table}: A data frame with ASV/OTUs as rows and samples 
#'       as columns.
#'     \item \code{metadata}: A data frame with samples as rows and grouping 
#'       variables as columns.
#'     \item \code{otu_ranked}: A data frame with ranked taxa containing:
#'       \itemize{
#'         \item \code{otu}: A column with taxa names.
#'         \item \code{rank}: A column with the rank for each taxon.
#'       }
#'     \item \code{elbow_core}: Character vector of OTU IDs identified as core
#'       using the elbow method.
#'     \item \code{increase_core}: Character vector of OTU IDs identified as core
#'       using the increase method.
#' }
#' 
#' @param core_set Which core set to plot: "elbow" (default) or "increase".
#' @param group_var Metadata column for bar coloring. Default: "sampling_date".
#' @param plot_type Allow to select between 3 different plot types: `bar` plot, `line` plot, or `heatmap`.
#'
#' @return A ggplot2 object that can be further customized.
#'
#' @examples
#' \dontrun{
#' # Generate an object from identify_core and then plot
#' data("switchgrass", package = "BRCore")
#' 
#' switchgrass_core <- identify_core(
#'   physeq_obj = switchgrass,
#'   priority_var = "sampling_date",  
#'   increase_value = 0.02,
#'   abundance_weight = 0,
#'   seed = 1234)
#' 
#' plot_core_distribution(core_result = switchgrass_core,
#'                        core_set = "increase",
#'                        group_var = "sampling_date",
#'                        plot_type  = "bar")
#' }
#'
#' @export
plot_core_distribution <- function(core_result,
                                   core_set   = "elbow",
                                   group_var  = "sampling_date",
                                   plot_type  = c("bar", "line", "heatmap")) {
    # match plot_type
    plot_type <- match.arg(plot_type)
    
    # Validate core_set argument
    if (!core_set %in% c("elbow", "increase")) {
        stop("core_set must be either 'elbow' or 'increase'")
    }
    
    # Select appropriate core OTU set
    if (core_set == "elbow") {
        core_otus <- core_result$elbow_core
    } else {
        core_otus <- core_result$increase_core
    }
    
    # Calculate relative abundance
    otu_relabun <- vegan::decostand(as.data.frame(core_result$otu_table),
                                    method = "total",
                                    MARGIN = 2)
    
    # Prepare metadata
    map <- as.data.frame(core_result$metadata)
    
    if (!"sample_id" %in% colnames(map)) {
        map <- map %>% tibble::rownames_to_column("sample_id")
    }
    
    # Validate required columns
    if (!group_var %in% colnames(map)) {
        stop(paste("Missing required column", group_var, "in metadata"))
    }
    
    otu_ranked <- core_result$otu_ranked
    
    # Prepare plotting data
    plotDF <- otu_relabun %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "otu") %>%
        tidyr::pivot_longer(cols = -otu,
                            names_to = "sample_id",
                            values_to = "relabun") %>%
        dplyr::left_join(map, by = "sample_id") %>%
        dplyr::left_join(otu_ranked, by = 'otu') %>%
        dplyr::filter(otu %in% core_otus) %>%
        dplyr::group_by(otu, .data[[group_var]]) %>%
        dplyr::summarise(
            time_freq = sum(relabun > 0) / dplyr::n(),
            coreTime = ifelse(time_freq == 1, 1, 0),
            detect = ifelse(time_freq > 0, 1, 0),
            .groups = "drop"
        )
    
    # Set factor levels for proper ordering
    #plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu)
    core_levels <- otu_ranked$otu[otu_ranked$otu %in% core_otus]
    plotDF$otu <- factor(plotDF$otu, levels = core_levels)
    
    # Create plot
    if (plot_type == "bar") {
        
        p <- ggplot2::ggplot(
            plotDF,
            ggplot2::aes(
                x    = otu,
                y    = time_freq,
                fill = factor(.data[[group_var]])
            )
        ) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::coord_flip() +
            ggplot2::scale_x_discrete(limits = rev(levels(plotDF$otu))) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                plot.title       = ggplot2::element_text(
                    hjust = 0.5, size = 12, face = "bold"
                ),
                plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 9),
                axis.text        = ggplot2::element_text(size = 6),
                legend.key.height = grid::unit(0.4, "cm"),
                legend.key.width  = grid::unit(0.4, "cm"),
                legend.title      = ggplot2::element_blank(),
                legend.text       = ggplot2::element_text(size = 8)
            ) +
            ggplot2::labs(
                title = paste("Core set occupancy across:", group_var),
                x     = "Ranked ASV/OTUs",
                y     = "Occupancy"
            ) +
            ggsci::scale_fill_npg()
        
    } else if (plot_type == "line") {
        
        group_levels <- length(levels(factor(plotDF[[group_var]])))
        
        # one facet per group_var level, one line per panel
        p <- ggplot2::ggplot(
            plotDF,
            ggplot2::aes(
                x     = otu,
                y     = time_freq,
                group = factor(.data[[group_var]]), 
                color = factor(.data[[group_var]])
            )
        ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 0.7) +
            ggplot2::facet_wrap(stats::as.formula(paste("~", group_var)), nrow = group_levels) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                plot.title        = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
                legend.position   = "none",
                axis.ticks.x      = element_line(),
                axis.text.y       = ggplot2::element_text(size = 8),
                axis.text.x       = ggplot2::element_text(angle = 45, size = 6, vjust = 1, hjust = 1)
            ) +
            ggplot2::labs(
                title = paste("Core set occupancy across:", group_var),
                x     = "Ranked ASV/OTUs",
                y     = "Occupancy"
            )
        
    } else if (plot_type == "heatmap") {
        
        p <- ggplot2::ggplot(
            plotDF,
            ggplot2::aes(
                x    = factor(.data[[group_var]]),
                y    = otu,
                fill = time_freq
            )
        ) +
            ggplot2::geom_tile(color = "white", size = 0.5) +
            # viridis if you want:
            # viridis::scale_fill_viridis(option = "plasma", name = "Occupancy") +
            ggplot2::theme_classic() +
            ggplot2::theme(
                plot.title        = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
                legend.key.height = grid::unit(0.4, "cm"),
                legend.key.width  = grid::unit(0.4, "cm"),
                legend.title      = ggplot2::element_text(size = 8, face = "bold"),
                legend.text       = ggplot2::element_text(size = 8),
                axis.ticks.x      = element_line(),
                axis.text.y       = ggplot2::element_text(size = 6),
                axis.text.x       = ggplot2::element_text(angle = 33, vjust = 1, hjust = 1)
            ) +
            ggplot2::labs(
                title = paste("Core set occupancy across:", group_var),
                x     = paste(group_var),
                y     = "Ranked ASV/OTUs",
                fill = "Occupancy"
            )
    }
    
    return(p)
}