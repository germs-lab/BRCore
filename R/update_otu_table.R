#' Add a rarefied otu_table to a phyloseq object
#'
#' This function imports a rarefied otu_table (or normalized) in a `phyloseq`
#' object `otu_table()` layer. It handles both summarized output (from
#' `multi_rarefy(.summarize = TRUE)`) and unsummarized output with iteration
#' suffixes (from `multi_rarefy(.summarize = FALSE)`).
#'
#' When unsummarized output is provided, the `sample_data` is replicated to
#' match the iteration rows. Each row retains its `_iter_#` suffix so
#' individual iterations can be distinguished downstream.
#'
#' @param physeq A `phyloseq` object in which you want to add the rarefied OTU/ASV table.
#' @param otu_rare A rarefied otu_table dataframe. Can be either:
#'   - Summarized: one row per sample (from `multi_rarefy(.summarize = TRUE)`)
#'   - Unsummarized: multiple rows per sample with `_iter_#` suffixes
#'     (from `multi_rarefy(.summarize = FALSE)`). The `sample_data` will be
#'     replicated to match each iteration row.
#'
#' @return A `phyloseq` object.
#'
#' @importFrom phyloseq otu_table sample_data tax_table phy_tree refseq
#' @importFrom phyloseq phyloseq sample_names prune_taxa prune_samples
#' @importFrom phyloseq sample_sums taxa_sums
#' @importFrom cli cli_inform cli_alert_success cli_alert_warning cli_alert_info
#' @importFrom magrittr %>%
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Perform multiple rarefaction (summarized)
#' otu_table_rare <-
#'   multi_rarefy(
#'     physeq = GlobalPatterns,
#'     depth_level = 200,
#'     num_iter = 3,
#'     threads = 1,
#'     set_seed = 123
#'   )
#'
#' rarefied_GlobalPatterns <-
#'   update_otu_table(
#'     physeq = GlobalPatterns,
#'     otu_rare = otu_table_rare
#'   )
#'
#' sample_sums(rarefied_GlobalPatterns)
#'
#' # Also works with unsummarized output (keeps all iterations)
#' otu_table_rare_full <-
#'   multi_rarefy(
#'     physeq = GlobalPatterns,
#'     depth_level = 200,
#'     num_iter = 3,
#'     .summarize = FALSE,
#'     threads = 1,
#'     set_seed = 123
#'   )
#'
#' rarefied_GlobalPatterns2 <-
#'   update_otu_table(
#'     physeq = GlobalPatterns,
#'     otu_rare = otu_table_rare_full
#'   )
#'
#' sample_sums(rarefied_GlobalPatterns2)
#' }
#'
#' @export
update_otu_table <- function(physeq, otu_rare) {
  physeq_samples <- sample_names(physeq)

  # Detect unsummarized output (row names contain _iter_ suffixes)
  otu_rare_rownames <- rownames(otu_rare)
  has_iter_suffix <- any(grepl("_iter_\\d+$", otu_rare_rownames))

  if (has_iter_suffix) {
    cli::cli_alert_info(
      "Detected unsummarized rarefaction output (iteration suffixes found). Keeping all iterations."
    )

    # Extract base sample names from iteration row names
    base_sample_ids <- sub("_iter_\\d+$", "", otu_rare_rownames)
    unique_samples <- unique(base_sample_ids)

    # Identify shared and removed samples using base names
    shared_samples <- intersect(physeq_samples, unique_samples)
    removed_samples <- setdiff(physeq_samples, unique_samples)

    .report_sample_status(
      shared_samples,
      removed_samples,
      otu_rare,
      physeq_samples
    )

    # Keep only rows whose base sample is in shared_samples
    keep_rows <- base_sample_ids %in% shared_samples
    otu_rare_ord <- otu_rare[keep_rows, , drop = FALSE]
    iter_rownames <- rownames(otu_rare_ord)
    iter_base_ids <- sub("_iter_\\d+$", "", iter_rownames)

    cli::cli_alert_info(
      "Building phyloseq object with {.val {nrow(otu_rare_ord)}} rows ({.val {length(shared_samples)}} samples x iterations) and {.val {ncol(otu_rare_ord)}} taxa"
    )

    # Replicate sample_data to match iteration rows
    sd_orig <- as(sample_data(physeq), "data.frame")
    sd_expanded <- sd_orig[iter_base_ids, , drop = FALSE]
    rownames(sd_expanded) <- iter_rownames

    # Build components
    phyloseq_components <- list(
      otu_table(
        t(otu_rare_ord) |> as.matrix() |> as.data.frame(),
        taxa_are_rows = TRUE
      ),
      sample_data(sd_expanded)
    )
  } else {
    # Summarized path
    otu_rare_samples <- rownames(otu_rare)
    shared_samples <- intersect(physeq_samples, otu_rare_samples)
    removed_samples <- setdiff(physeq_samples, otu_rare_samples)

    .report_sample_status(
      shared_samples,
      removed_samples,
      otu_rare,
      physeq_samples
    )

    otu_rare_ord <- otu_rare[shared_samples, , drop = FALSE]

    cli::cli_alert_info(
      "Building phyloseq object with {.val {nrow(otu_rare_ord)}} samples and {.val {ncol(otu_rare_ord)}} taxa"
    )

    phyloseq_components <- list(
      otu_table(
        t(otu_rare_ord) |> as.matrix() |> as.data.frame(),
        taxa_are_rows = TRUE
      ),
      sample_data(physeq)[shared_samples, ]
    )
  }

  # Add optional components if they exist
  .add_optional_component <- function(components, accessor) {
    obj <- accessor(physeq, errorIfNULL = FALSE)
    if (!is.null(obj)) c(components, list(obj)) else components
  }

  phyloseq_components <- .add_optional_component(phyloseq_components, tax_table)
  phyloseq_components <- .add_optional_component(phyloseq_components, phy_tree)
  phyloseq_components <- .add_optional_component(phyloseq_components, refseq)

  # Build the new phyloseq object
  new_phyloseq <- do.call(phyloseq, phyloseq_components) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    prune_samples(sample_sums(x = .) > 0, x = .)

  cli::cli_alert_success("Analysis complete!")

  return(new_phyloseq)
}


#' Report sample status after rarefaction (internal helper)
#' @noRd
.report_sample_status <- function(
  shared_samples,
  removed_samples,
  otu_rare,
  physeq_samples
) {
  if (length(shared_samples) == 0) {
    cli::cli_alert_danger(
      "No matching sample names between phyloseq object and rarefied OTU table."
    )
    stop("No matching sample names. Check rownames of otu_rare.")
  }

  otu_rare_samples <- rownames(otu_rare)

  if (identical(physeq_samples, otu_rare_samples)) {
    cli::cli_alert_success(
      "Phyloseq object and rarefied otu_tables sample names are identical."
    )
  } else {
    cli::cli_alert_warning(
      "Phyloseq object and rarefied otu_tables sample names are NOT identical. Check below samples removed by rarefaction."
    )
  }

  if (length(removed_samples) > 0) {
    cli::cli_alert_warning(
      "{.val {length(removed_samples)}} sample{?s} removed due to rarefaction: {.val {paste(removed_samples, collapse = ', ')}}"
    )
  } else {
    rarefaction_depth <- sum(otu_rare[1, ])
    cli::cli_alert_success(
      "All samples kept after rarefaction at depth of: {.val {rarefaction_depth}}"
    )
  }
}
