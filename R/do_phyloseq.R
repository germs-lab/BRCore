#' Add a rarefied otu_table to a phyloseq object
#'
#' This function imports a rarefied otu_table (or normalized) in a `phyloseq` object `otu_table()` layer.
#'
#' @param physeq A `phyloseq` object in which you want to add the rarefied OTU/ASV table.
#' @param otu_rare A rarefied otu_table dataframe.
#'
#' @return A `phyloseq` object.
#'
#' @import phyloseq
#' @import vegan
#' @import dplyr
#' @import tibble
#' @importFrom dplyr %>% 
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Perform multiple rarefaction
#' otu_table_rare <- 
#'     multi_rarefy(
#'         physeq = GlobalPatterns,
#'         depth_level = 200,
#'         num_iter = 3,
#'         threads = 1,
#'         set_seed = 123
#'     )
#'
#' # Check the rarefied data output
#' rowSums(otu_table_rare)
#'
#' # Recreate the phyloseq object and check
#' rarefied_GlobalPatterns<- 
#'     do_phyloseq(physeq = GlobalPatterns, 
#'                 otu_rare = otu_table_rare )
#'     
#' sample_sums(rarefied_GlobalPatterns)
#'}
#'
#' @export
do_phyloseq <- function(physeq, 
                        otu_rare) {
    
    # sample names from original phyloseq object
    physeq_samples <- sample_names(physeq)
    cli::cli_inform("Sample names from the original phyloseq object (otu_table):")
    print(physeq_samples)
    
    # reorder rarefied OTU table to match
    otu_rare_ord <- otu_rare[physeq_samples, , drop = FALSE]
    
    # sample names from rarefied OTU table (after reordering)
    otu_rare_samples <- rownames(otu_rare_ord)
    cli::cli_inform("Sample names in the rarefied OTU table (data.frame) after reordering to match phyloseq:")
    print(otu_rare_samples)
    
    if (identical(physeq_samples, otu_rare_samples)) {
        cli::cli_alert_success("Phyloseq object and rarefied otu_tables sample names are identical.")
    } else {
        cli::cli_alert_warning("Phyloseq object and rarefied otu_tables sample names are NOT identical. Check below samples removed by rarefaction.")
    }
    
    # Check for samples removed due to rarefaction
    if (!all(physeq_samples %in% otu_rare_samples)) {
        removed_samples <- physeq_samples[!(physeq_samples %in% otu_rare_samples)]
        cli::cli_alert_warning("Samples removed due to rarefaction: {paste(removed_samples, collapse = ', ')}")
    } else {
        # Calculate rarefaction depth (assuming all samples have the same depth after rarefaction)
        rarefaction_depth <- sum(otu_rare[1, ])
        cli::cli_alert_success("All samples kept after rarefaction at depth of: {rarefaction_depth}")
    }
    
    # Create list of components for the new phyloseq object
    phyloseq_components <- list(
        phyloseq::otu_table(
            t(otu_rare_ord) %>%
                as.matrix() %>%
                as.data.frame(),
            taxa_are_rows = TRUE
        ),
        phyloseq::sample_data(physeq)
    )
    
    # Add optional components if they exist and are not empty
    if (!is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
        phyloseq_components <- c(phyloseq_components, list(phyloseq::tax_table(physeq)))
    }
    
    if (!is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))) {
        phyloseq_components <- c(phyloseq_components, list(phyloseq::phy_tree(physeq)))
    }
    
    if (!is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))) {
        phyloseq_components <- c(phyloseq_components, list(phyloseq::refseq(physeq)))
    }
    
    # Build the new phyloseq object
    new_phyloseq <- do.call(phyloseq::phyloseq, phyloseq_components) %>%
        phyloseq::prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        phyloseq::prune_samples(sample_sums(x = .) > 0, x = .)
    
    cli::cli_alert_success("Analysis complete!")
    
    return(new_phyloseq)
}