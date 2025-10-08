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
#'                        group_var = "sampling_date")
#' 
#' }
#'
#' @export
plot_core_distribution <- function(core_result,
                                   core_set = "elbow",
                                   group_var = "sampling_date") {
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
    p <- ggplot(plotDF, aes(
        x = otu,
        y = time_freq,
        fill = factor(.data[[group_var]])
    )) +
        geom_bar(stat = 'identity', position = 'dodge') +
        coord_flip() +
        scale_x_discrete(limits = rev(levels(plotDF$otu))) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5,
                size = 12,
                face = "bold"
            ),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
            axis.text = element_text(size = 6),
            legend.key.height = grid::unit(0.4, "cm"),
            legend.key.width = grid::unit(0.4, "cm"),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 8)
        ) +
        labs(
            title = paste("Core set occupancy across:", group_var),
            x = "Ranked ASV/OTUs",
            y = "Occupancy"
        ) +
        ggsci::scale_fill_npg()
    
    return(p)
}
