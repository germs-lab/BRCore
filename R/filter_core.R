#' Select ASVs Based on Occurrence Threshold
#'
#' This function selects ASVs from a phyloseq object based on their occurrence frequency
#' across samples. It returns filtered phyloseq objects for high and low occurrence ASVs.
#'
#' @param physeq A phyloseq object.
#' @param threshold A numeric value between 0 and 1 specifying the occurrence threshold.
#'   Default is 0.6 (60% of samples).
#' @param as The matrix dimension that is desired for taxa. Must be
#'   "rows" for rows and "columns" or "cols" for columns.
#' @return A named list containing:
#'   \itemize{
#'     \item \code{physeq_high_occ}: Phyloseq object with high occurrence ASVs.
#'     \item \code{physeq_low_occ}: Phyloseq object with low occurrence ASVs.
#'   }
#' @import phyloseq
#' @import dplyr
#' @export
filter_core <- function(physeq, threshold = 0.6, as = "rows") {
  # Validate inputs
  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort("{.arg physeq} must be a phyloseq object. Got {class(physeq)[1]} instead.")
  }

  if (threshold < 0 || threshold > 1) {
    cli::cli_abort("{.arg threshold} must be between 0 and 1. Got {threshold} instead.")
  }

  # Convert to relative abundance
  # physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

  # Helper function to select taxa orientation
  # Based on `orient_taxa()` https://github.com/mikemc/speedyseq/blob/0057652ff7a4244ccef2b786dca58d901ec2fc62/R/utils.R
  ensure_phyloseq_orientation <- function(physeq, as) {
    # Validate 'as' argument
    stopifnot(
      "Arguments 'as' must be a string:'rows', 'cols' or 'columns'!" =
        as %in% c("rows", "columns", "cols")
    )

    # Ensure taxa are in the desired orientation
    if (identical(as, "rows")) {
      if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
      }
    } else {
      if (taxa_are_rows(physeq)) {
        physeq <- t(physeq)
      }
    }

    return(physeq)
  }
  physeq_rel <- ensure_phyloseq_orientation(physeq, as)

  # Calculate occurrence of each ASV across samples
  # Get total number of samples
  if (identical(as, "rows")) {
    asv_sample_counts <- rowSums(otu_table(physeq_rel) > 0)
    total_samples <- ncol(otu_table(physeq_rel))
  }

  if (identical(as, "cols") || identical(as, "columns")) {
    asv_sample_counts <- colSums(otu_table(physeq_rel) > 0)
    total_samples <- nrow(otu_table(physeq_rel))
  }

  # Select ASVs based on occurrence criteria
  high_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts >= total_samples * threshold])
  low_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts < total_samples * threshold])

  # Print ASV counts
  cli::cli_alert_info("ASVs found in >={threshold * 100}% samples: {length(high_occurrence_asvs)}")
  cli::cli_alert_info("ASVs found in <{threshold * 100}% samples: {length(low_occurrence_asvs)}")

  # Filter phyloseq objects for each category
  physeq_high_occ <- prune_taxa(high_occurrence_asvs, physeq_rel)
  physeq_low_occ <- prune_taxa(low_occurrence_asvs, physeq_rel)

  # Return results
  return(list(
    physeq_high_occ = physeq_high_occ,
    physeq_low_occ = physeq_low_occ
  ))
}
