#' Subset phyloseq object based on core analysis results
#'
#' Filters a phyloseq object to retain only samples and taxa identified in a core microbiome analysis
#' using `extract_core()`.
#' Returns both the subsetted phyloseq object and corresponding ASV matrix.
#'
#' @param core_obj A list object containing core microbiome analysis results.
#'        Expected to contain a dataframe at position \code{[[4]]} with OTU information.
#' @param physeq A phyloseq object to be subsetted.
#' @param .var Character string specifying the column name in \code{core_obj[[4]]}
#'        containing sample/group identifiers to filter by.
#' @param type Value to filter by in the 'fill' column of \code{core_obj[[4]]}.
#'        Determines which OTUs are retained.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{expected_physeq} - Subsetted phyloseq object
#'   \item \code{asv_matrix} - Corresponding ASV abundance matrix
#' }
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#'   \item Extracts OTUs of interest from the core analysis object based on the specified type
#'   \item Filters the phyloseq object to retain only these OTUs and samples where they appear
#'   \item Verifies the filtering was successful
#'   \item Returns both the subsetted phyloseq object and ASV matrix
#' }
#'
#' The function uses \code{prune_samples} and \code{prune_taxa} to safely subset the phyloseq object
#' while preserving all associated data (sample data, taxonomy table, etc.). The ASV matrix is
#' processed to remove samples with zero row sums.
#'
#' @examples
#' \dontrun{
#' # Subset for "persistent" core OTUs
#' result <- subset_physeq(
#'   core_obj = my_core_results,
#'   physeq = my_phyloseq,
#'   .var = "body_site",
#'   type = "persistent"
#' )
#'
#' # Access results
#' subset_physeq <- result$expected_physeq
#' asv_matrix <- result$asv_matrix
#' }
#'
#' @import phyloseq
#' @import dplyr
#' @import tibble
#' @import cli
#' @export
subset_physeq <- function(core_obj, physeq, .var, type) {
  asv_matrix <-
    core_obj[[4]] %>% # Get strings from all OTUs
    dplyr::filter(., .$fill == type) %>%
    tibble::column_to_rownames(., var = .var) %>%
    rownames() %>%
    extract_matrix(physeq, .vec = .) # Remove samples where row sum 0 through extract_matrix()

  sample_strings <- rownames(asv_matrix)

  ## Phyloseqs
  expected_physeq <- phyloseq::prune_samples(
    sort(phyloseq::sample_names(physeq)) %in% sort(sample_strings),
    physeq
  ) %>%
    phyloseq::prune_taxa(rownames(.@otu_table) %in% colnames(asv_matrix), .)

  ## Checking that core is filtered out
  if (any(rownames(expected_physeq@otu_table) %in% colnames(asv_matrix))) {
    cli::cli_alert_info("Successfully selected ASVs of interest")
  }

  return(
    list(
      expected_physeq = expected_physeq,
      asv_matrix = asv_matrix
    )
  )
}
