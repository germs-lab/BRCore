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
#' @examplesIf requireNamespace("phyloseq", quietly = TRUE)
#' \donttest{
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
#' @importFrom magrittr %>%
#' @importFrom grid unit
#' @export
plot_identified_core <- function(bray_curtis_ranked, 
                                 elbow, 
                                 lastCall,
                                 increase_value = 0.02) {
    
    new_data <- bray_curtis_ranked %>%
        dplyr::mutate(rank_num = suppressWarnings(as.integer(as.character(rank)))) %>%
        dplyr::filter(rank_num <= 1.2 * lastCall) %>%
        dplyr::mutate(rank_fac = factor(rank_num, levels = rank_num))
    
    # Calculate percentage label text
    percent_label <- paste0(round(increase_value * 100), "%")
    
    # Calculate dynamic label positions
    max_rank <- max(new_data$rank_num, na.rm = TRUE)
    max_prop <- max(new_data$proportionBC, na.rm = TRUE)
    min_prop <- min(new_data$proportionBC, na.rm = TRUE)
    prop_range <- max_prop - min_prop
    
    # Dynamic y positions for labels (avoid overlap)
    elbow_y <- min_prop + 0.15 * prop_range
    lastCall_y <- min_prop + 0.35 * prop_range
    
    # Adjust x position for labels if too close to edges
    elbow_x_offset <- ifelse(elbow < max_rank * 0.1, elbow + max_rank * 0.02, elbow - max_rank * 0.02)
    lastCall_x_offset <- ifelse(lastCall > max_rank * 0.9, lastCall - max_rank * 0.02, lastCall + max_rank * 0.02)
    
    # Create base plot
    p <- ggplot2::ggplot(new_data, ggplot2::aes(x = rank_num, y = proportionBC)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.7) +
        
        # Enhanced vertical lines with better visibility
        ggplot2::geom_vline(xintercept = elbow, linetype = "dashed", 
                            color = "#d73027", linewidth = 1.2, alpha = 0.8) +
        ggplot2::geom_vline(xintercept = lastCall, linetype = "dashed", 
                            color = "#1a9850", linewidth = 1.2, alpha = 0.8) +
        
        # Add colored background rectangles to highlight core regions
        ggplot2::annotate("rect", xmin = 0, xmax = elbow, 
                          ymin = -Inf, ymax = Inf, 
                          fill = "white", alpha = 0.05) +
        ggplot2::annotate("rect", xmin = elbow, xmax = lastCall, 
                          ymin = -Inf, ymax = Inf, 
                          fill = "white", alpha = 0.05) +
        
        # Improved text labels with backgrounds and better positioning
        ggplot2::annotate("label", x = elbow_x_offset, y = elbow_y, 
                          label = paste0("Elbow (", elbow, ")"), 
                          color = "#d73027", fill = "white", 
                          size = 3.5, fontface = "bold",
                          label.padding = ggplot2::unit(0.3, "lines")) +
        ggplot2::annotate("label", x = lastCall_x_offset, y = lastCall_y, 
                          label = paste0("Last ", percent_label, " (", lastCall, ")"), 
                          color = "#1a9850", fill = "white", 
                          size = 3.5, fontface = "bold",
                          label.padding = ggplot2::unit(0.3, "lines")) +
        
        # Enhanced styling
        ggplot2::labs(
            title = paste0("Core Microbiome Identification: Elbow and Last ", percent_label, " Methods"),
            subtitle = paste0("Core taxa identified: Elbow method (", elbow, " OTUs), Last ", 
                              percent_label, " method (", lastCall, " OTUs)"),
            x = "Ranked OTUs", 
            y = "% Bray-Curtis similarity"
        )  +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)) +
      # Ensure proper scaling
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05)))
    
    # To add some more fancy graphics later...
    # axis.title = ggplot2::element_text(size = 11),
    # axis.text = ggplot2::element_text(size = 9),
    # panel.grid.minor = ggplot2::element_blank(),
    # panel.grid.major.x = ggplot2::element_line(color = "gray90", size = 0.3),
    # panel.grid.major.y = ggplot2::element_line(color = "gray90", size = 0.3))
    
    return(p)
}
