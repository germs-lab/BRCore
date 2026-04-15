#' Add a rarefied otu_table to a phyloseq object
#'
#' This function updates a `phyloseq` object by replacing its OTU/ASV table
#' with a rarefied version produced by `multi_rarefy()`. The rarefied table can
#' be a data frame, a list of data frames (`.as_array = FALSE`), or a 3D array
#' (`.as_array = TRUE`). When providing a list or array, specify which iteration
#' to use via the `iteration` parameter.
#'
#' @param physeq_obj A `phyloseq` object in which you want to add the rarefied
#' OTU/ASV table.
#' @param rarefied_otus A data frame, list of data frames, or 3D array output from
#' `multi_rarefy()` containing the rarefied OTU/ASV tables.
#' @param iteration Integer specifying which iteration to extract from a list
#' or array. Required when `rarefied_otus` is a list or array. Ignored when
#' `rarefied_otus` is a data frame.
#'
#' @return A `phyloseq` object.
#'
#' @importFrom phyloseq otu_table sample_data tax_table phy_tree refseq
#' @importFrom phyloseq phyloseq sample_names prune_taxa prune_samples
#' @importFrom phyloseq sample_sums taxa_sums
#' @importFrom cli cli_inform cli_alert_success cli_alert_warning cli_alert_info
#' @importFrom cli cli_abort
#' @importFrom magrittr %>%
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # List output (.as_array = FALSE)
#' otu_list <-
#'   multi_rarefy(
#'     physeq_obj = GlobalPatterns,
#'     depth_level = 200,
#'     num_iter = 3,
#'     .as_array = FALSE,
#'     set_seed = 123
#'   )
#'
#' # Extract iteration 2
#' rarefied_gp <- update_otu_table(GlobalPatterns, otu_list, iteration = 2)
#'
#' # Array output (.as_array = TRUE)
#' otu_array <-
#'   multi_rarefy(
#'     physeq_obj = GlobalPatterns,
#'     depth_level = 200,
#'     num_iter = 3,
#'     .as_array = TRUE,
#'     set_seed = 123
#'   )
#'
#' # Extract iteration 1
#' rarefied_gp2 <- update_otu_table(GlobalPatterns, otu_array, iteration = 1)
#' }
#'
#' @export
update_otu_table <- function(physeq_obj, rarefied_otus, iteration = NULL) {
  suppressMessages(.phyloseq_class_check(physeq_obj))

  physeq_samples <- sample_names(physeq_obj)

  # Handle list input
  if (is.list(rarefied_otus) && !is.data.frame(rarefied_otus)) {
    if (is.null(iteration)) {
      cli::cli_abort(
        "{.arg iteration} must be specified when {.arg rarefied_otus} is a list."
      )
    }

    n_iter <- length(rarefied_otus)

    if (iteration < 1 || iteration > n_iter) {
      cli::cli_abort(
        "{.arg iteration} must be between 1 and {n_iter}."
      )
    }

    cli::cli_alert_info(
      "Extracting iteration {.val {iteration}} from list of {.val {n_iter}} iterations."
    )

    rarefied_otus <- as.data.frame(rarefied_otus[[iteration]])
  }

  # Handle 3D array input
  if (is.array(rarefied_otus) && length(dim(rarefied_otus)) == 3) {
    if (is.null(iteration)) {
      cli::cli_abort(
        "{.arg iteration} must be specified when {.arg rarefied_otus} is a 3D array."
      )
    }

    n_iter <- dim(rarefied_otus)[3]

    if (iteration < 1 || iteration > n_iter) {
      cli::cli_abort(
        "{.arg iteration} must be between 1 and {n_iter}."
      )
    }

    cli::cli_alert_info(
      "Extracting iteration {.val {iteration}} from array with {.val {n_iter}} iterations."
    )

    rarefied_otus <- as.data.frame(rarefied_otus[,, iteration])
  }

  # Now rarefied_otus is always a data frame
  otu_rare_samples <- rownames(rarefied_otus)
  shared_samples <- intersect(physeq_samples, otu_rare_samples)
  removed_samples <- setdiff(physeq_samples, otu_rare_samples)

  .report_sample_status(
    shared_samples,
    removed_samples,
    rarefied_otus,
    physeq_samples
  )

  otu_rare_ord <- rarefied_otus[shared_samples, , drop = FALSE]

  cli::cli_alert_info(
    "Building phyloseq object with {.val {nrow(otu_rare_ord)}} samples and {.val {ncol(otu_rare_ord)}} taxa"
  )

  # Build components
  phyloseq_components <- list(
    otu_table(
      t(otu_rare_ord) |> as.matrix() |> as.data.frame(),
      taxa_are_rows = TRUE
    ),
    sample_data(physeq_obj)[shared_samples, ]
  )

  # Add optional components if they exist
  .add_optional_component <- function(components, accessor) {
    obj <- accessor(physeq_obj, errorIfNULL = FALSE)
    if (!is.null(obj)) c(components, list(obj)) else components
  }

  phyloseq_components <- .add_optional_component(phyloseq_components, tax_table)
  phyloseq_components <- .add_optional_component(phyloseq_components, phy_tree)
  phyloseq_components <- .add_optional_component(phyloseq_components, refseq)

  # Build the new phyloseq object
  new_phyloseq <- do.call(phyloseq, phyloseq_components) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    prune_samples(sample_sums(x = .) > 0, x = .)

  cli::cli_alert_success("Update complete!")

  return(new_phyloseq)
}


#' Report sample status after rarefaction
#' Internal function to report on the status of samples after rarefaction, including any removed samples and whether sample names match between the original phyloseq object and the rarefied OTU table.
#' @noRd
#' @keywords internal
.report_sample_status <- function(
  shared_samples,
  removed_samples,
  rarefied_otus,
  physeq_samples
) {
  if (length(shared_samples) == 0) {
    cli::cli_alert_danger(
      "No matching sample names between phyloseq object and rarefied OTU table."
    )
    stop("No matching sample names. Check rownames of rarefied_otus.")
  }

  otu_rare_samples <- rownames(rarefied_otus)

  if (identical(sort(physeq_samples), sort(otu_rare_samples))) {
    cli::cli_alert_success(
      "Phyloseq object and rarefied otu_table sample names are identical."
    )
  } else {
    cli::cli_alert_warning(
      "Phyloseq object and rarefied otu_table sample names are NOT identical. Check samples removed by rarefaction below."
    )
  }

  if (length(removed_samples) > 0) {
    cli::cli_alert_warning(
      "{.val {length(removed_samples)}} sample{?s} removed due to rarefaction: {.val {paste(removed_samples, collapse = ', ')}}"
    )
  } else {
    rarefaction_depth <- sum(rarefied_otus[1, ])
    cli::cli_alert_success(
      "All samples kept after rarefaction at depth of: {.val {rarefaction_depth}}"
    )
  }
}
