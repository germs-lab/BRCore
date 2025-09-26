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
#'     parallel_rarefy(
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
            physeq@tax_table,
            physeq@refseq
        ) %>%
        prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        prune_samples(sample_sums(x = .) > 0, x = .)

    return(new_phyloseq)
}
