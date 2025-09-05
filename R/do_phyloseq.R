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
#' # Load example data
#' library(phyloseq)
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Perform multiple rarefaction
#' rarefied_data <- 
#'     multi_rarefy(GlobalPatterns, 
#'     depth_level = 500, 
#'     num_iter = 10, 
#'     threads = 2, 
#'     set_seed = 453)
#'
#' # Check the rarefied data output
#' rowSums(rarefied_data)
#'
#' if (is.null(rarefied_data) || nrow(rarefied_data) == 0 || ncol(rarefied_data) == 0) {
#' stop("rarefied_data has zero dimensions. Check data preprocessing and input files.")
#' }
#'
#' # Recreate the phyloseq object and check
#' rarefied_physeq <- 
#'     do_phyloseq(physeq = GlobalPatterns, 
#'                 otu_rare = rarefied_data )
#'     
#' rarefied_physeq
#'
#' @export
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
