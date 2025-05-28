#' Import a rarefied otu_table into a phyloseq object
#'
#' This function imports a rarefied otu_table (or normalized) in a `phyloseq` object `otu_table()` layer.
#'
#' @param physeq A `phyloseq` object in which I want to add the rarefied OTU/ASV table.
#' @param otu_rare A rarefied otu_table dataframe.
#'
#' @return A `phyloseq` object. 
#' 
#' @import phyloseq
#' @import vegan
#' @import dplyr
#' @import tibble
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Perform multiple rarefaction
#' rarefied_data <- multi_rarefy(GlobalPatterns, depth_level = 1000, num_iter = 99)
#' head(rarefied_data)
#' 
#' rarefied_physeq <- do_phyloseq(physeq = GlobalPatterns, otu_rare=rarefied_data )
#' head(otu_table(rarefied_physeq))
#' }
do_phyloseq <- function(physeq, otu_rare) {
    new_phyloseq <-
        phyloseq(
            otu_table(
                otu_rare %>%
                    t() %>%
                    as.matrix() %>%
                    as.data.frame(),
                taxa_are_rows = TRUE
            ),
            physeq@sam_data,
            physeq@tax_table
        ) %>%
        prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        prune_samples(sample_sums(x = .) > 0, x = .)
    
    return(new_phyloseq)
}
