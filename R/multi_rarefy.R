#' Run multiple rarefaction for microbiome count tables
#'
#' This function performs multiple rarefaction on a `phyloseq` object by randomly
#' sub-sampling OTUs/ASVs within samples without replacement. The process is
#' repeated for a specified number of iterations, and the results are averaged.
#' Samples with fewer OTUs/ASVs than the specified `depth_level` are discarded.
#'
#' @param physeq_obj A `phyloseq` object containing an OTU/ASV table.
#' @param depth_level An integer specifying the sequencing depth (number of
#'   OTUs/ASVs) to which samples should be rarefied.
#' @param num_iter An integer specifying the number of iterations to perform
#'   for rarefaction.
#' @param .as A character string indicating whether to return the results as a
#' 3D array or as a list of data frames. If `"array"`, returns a 3D array with
#' dimensions (samples x taxa x iterations). If `"list"`, returns a list of
#' data frames, one for each iteration, with samples as rows and taxa as
#' columns. (default = "list")
#' @param set_seed An optional integer to set the random seed for
#' reproducibility (default = NULL).
#'
#' @return A data frame with taxa as rows and samples as columns. The values
#'   represent the average sequence counts calculated across all iterations.
#'   Samples with less than `depth_level` sequences are discarded.
#'
#'
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom parallel clusterEvalQ
#' @importFrom dplyr group_by summarise across everything filter near
#' @importFrom dplyr where
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table
#' @importFrom cli cli_h1 cli_h2 cli_alert_info cli_alert_warning cli_alert_success cli_alert_danger
#' @importFrom utils head
#' @importFrom vegan rrarefy
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' data("bcse", package = "BRCore")
#'
#' # Example rarefaction (single iteration, single core to keep examples fast)
#' otu_table_rare <- multi_rarefy(
#'   physeq_obj = bcse,
#'   depth_level = 1000,
#'   num_iter = 100,
#'   .as_array = FALSE,
#'   set_seed = 7642
#' )
#'
#' rowSums(otu_table_rare[[1]])
#' }
#'
#' @export
multi_rarefy <- function(
  physeq_obj,
  depth_level,
  num_iter = 100,
  .as = "list",
  set_seed = NULL
) {
  # Input validation ----
  cli::cli_h1("Multiple Rarefaction")
  cli::cli_h2("Input Validation")

  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    cli::cli_alert_danger(
      "The 'phyloseq' package is required but not installed."
    )
    stop("Please install 'phyloseq' to use this function.")
  }

  .phyloseq_class_check(physeq_obj)

  if (!is.numeric(depth_level) || depth_level <= 0) {
    cli::cli_alert_danger("depth_level must be a positive number")
    stop("Invalid depth_level")
  }

  if (!is.numeric(num_iter) || num_iter <= 0) {
    cli::cli_alert_danger("num_iter must be a positive integer")
    stop("Invalid num_iter")
  }

  # Seed setup: iteration-specific seeds deterministically ---
  if (is.null(set_seed)) {
    cli::cli_alert_warning("No seed set. Results may not be reproducible.")
    iteration_seeds <- rep(NA_integer_, num_iter)
  } else {
    cli::cli_alert_info("Seed: {.val {set_seed}}")
    set.seed(set_seed)
    iteration_seeds <- sample.int(.Machine$integer.max, num_iter)
  }

  # Prepare data ----
  otu_mat <- .extract_otu_matrix(physeq_obj, samples_as_rows = TRUE)

  ## Filter samples by depth ----
  n_samples_before <- nrow(otu_mat)
  original_sample_ids <- rownames(otu_mat) |>
    unique()
  n_taxa_before <- ncol(otu_mat)

  keep <- rowSums(otu_mat) >= depth_level

  if (length(keep) != nrow(otu_mat)) {
    stop("Logical subsetting mismatch: check OTU table orientation.")
  }

  otu_mat <- otu_mat[keep, , drop = FALSE]

  n_samp <- nrow(otu_mat)
  n_taxa <- ncol(otu_mat)

  # Input parameter checks ----
  cli::cli_alert_info(
    "Input (matrix/df dim): {.val {nrow(otu_mat)}} samples x {.val {ncol(otu_mat)}} taxa"
  )
  cli::cli_alert_info("Rarefaction depth: {.val {depth_level}}")
  cli::cli_alert_info("Iterations: {.val {num_iter}}")

  cli::cli_alert_info("taxa_are_rows: {phyloseq::taxa_are_rows(physeq_obj)}")

  cli::cli_alert_info(
    "OTU matrix/df rownames head: {paste(head(rownames(otu_mat)), collapse = ', ')}"
  )
  cli::cli_alert_info(
    "OTU matrix/df colnames head: {paste(head(colnames(otu_mat)), collapse = ', ')}"
  )

  cli::cli_alert_info(
    "Row sums summary: Min={min(rowSums(otu_mat))}, Max={max(rowSums(otu_mat))}, Median={median(rowSums(otu_mat))}"
  )

  # Single iteration ----
  if (num_iter == 1) {
    rare <- vegan::rrarefy(otu_mat, sample = depth_level)
    return(as.data.frame(rare))
  }

  # Multiple iterations ----
  if (.as == "array") {
    processed_data <- array(
      0,
      dim = c(n_samp, n_taxa, num_iter),
      dimnames = list(
        rownames(otu_mat),
        colnames(otu_mat),
        paste0("iter_", seq_len(num_iter))
      )
    )

    for (i in seq_len(num_iter)) {
      processed_data[,, i] <- vegan::rrarefy(otu_mat, sample = depth_level)
    }
  }
  if (.as == "list") {
    processed_data <- lapply(seq_len(num_iter), function(i) {
      if (!is.na(iteration_seeds[i])) {
        set.seed(iteration_seeds[i])
      }

      rare_result <- vegan::rrarefy(otu_mat, sample = depth_level)

      if (is.null(dim(rare_result))) {
        rare_result <- matrix(
          rare_result,
          nrow = 1,
          dimnames = list(rownames(otu_mat), names(rare_result))
        )
      }

      as.data.frame(rare_result)
    })
    names(processed_data) <- paste0("iter_", seq_len(num_iter))
  }

  avg_taxa_removed <- NULL

  if (.as == "array") {
    n_samples_removed <- n_samples_before - n_samp
    removed_samples <- setdiff(original_sample_ids, rownames(otu_mat))
    n_taxa_removed <- n_taxa - dim(processed_data)[2]
  }
  if (.as == "list") {
    n_taxa_after <- min(vapply(processed_data, ncol, integer(1)))

    n_samples_after <- nrow(processed_data[[1]])
    n_samples_removed <- n_samples_before - n_samples_after
    removed_samples <- setdiff(
      rownames(otu_mat),
      rownames(processed_data[[1]])
    )

    removed_samples <- setdiff(original_sample_ids, rownames(otu_mat))

    n_taxa_removed <- n_taxa_before - n_taxa_after

    avg_taxa_removed <- mean(vapply(processed_data, ncol, double(1)))
  }
  .report_rarefaction_results(
    n_samples_removed = n_samples_removed,
    n_taxa_before = n_taxa_before,
    n_taxa_removed = n_taxa_removed,
    avg_taxa_removed = avg_taxa_removed,
    removed_samples = removed_samples,
    dataframe = otu_mat,
    processed_data = processed_data,
    depth_level = depth_level,
    unique_samples = original_sample_ids,
    num_iter = num_iter
  )

  # Return the processed data frame
  return(processed_data)
}


#' Count zeros and sparsity in a data frame
#' @param dataframe A data frame to analyze.
#' @return A list containing the number of zeros
#'
#' @noRd
#' @keywords internal
.sparsity_count <- function(dataframe) {
  if (inherits(dataframe, "list")) {
    min_zeros <- min(vapply(dataframe, function(df) sum(df == 0), integer(1)))
    min_total <- min(vapply(
      dataframe,
      function(df) nrow(df) * ncol(df),
      integer(1)
    ))

    max_zeros <- max(vapply(dataframe, function(df) sum(df == 0), integer(1)))
    max_total <- max(vapply(
      dataframe,
      function(df) nrow(df) * ncol(df),
      integer(1)
    ))

    avg_zeros <- mean(vapply(dataframe, function(df) sum(df == 0), double(1)))
    avg_total <- mean(vapply(
      dataframe,
      function(df) nrow(df) * ncol(df),
      double(1)
    ))
  } else {
    original_zeros <- sum(dataframe == 0)
    original_total <- nrow(dataframe) * ncol(dataframe)
  }

  if (inherits(dataframe, "list")) {
    sparsity <- list(
      min_zeros = min_zeros,
      max_zeros = max_zeros,
      avg_zeros = avg_zeros,
      min_total = min_total,
      max_total = max_total,
      avg_total = avg_total,
      min_sparsity = round(min_zeros / min_total * 100, 2),
      max_sparsity = round(max_zeros / max_total * 100, 2),
      avg_sparsity = round(avg_zeros / avg_total * 100, 2)
    )
  } else {
    sparsity <- round(original_zeros / original_total * 100, 2)

    list(
      zeros = original_zeros,
      total = original_total,
      sparsity = sparsity
    )
  }
}

#' Helper function to report rarefaction results
#' @param n_samples_removed Number of samples removed during rarefaction.
#' @param n_taxa_removed Number of taxa removed during rarefaction.
#' @param removed_samples A vector of sample names that were removed.
#' @param dataframe The original data frame before rarefaction.
#' @param processed_data The data frame after rarefaction.
#' @param depth_level The depth level used for rarefaction.
#' @return None. This function is used for side-effect reporting via cli messages.
#' @noRd
#' @keywords internal
#'
.report_rarefaction_results <- function(
  n_samples_removed,
  n_taxa_before,
  n_taxa_removed,
  avg_taxa_removed,
  removed_samples,
  dataframe,
  processed_data,
  depth_level,
  unique_samples = NULL,
  num_iter = NULL
) {
  cli::cli_h2("Rarefaction Results")

  # Detect if processed_data is a 3D array or list of data frames

  is_array <- is.array(processed_data) && length(dim(processed_data)) == 3

  # Sample Removal
  cli::cli_h3("Sample Removal")
  if (n_samples_removed > 0) {
    cli::cli_alert_warning(
      "{.val {n_samples_removed}} sample{?s} removed (depth < {.val {depth_level}})"
    )
    cli::cli_alert_warning(
      "Samples removed: {.val {paste(removed_samples, collapse = ', ')}}"
    )
  } else {
    cli::cli_alert_success("No samples removed.")
  }

  # Taxa Removal
  cli::cli_h3("Taxa Removal")
  if (n_taxa_removed == 0) {
    cli::cli_alert_success("No taxa removed.")
    cli::cli_alert_warning(
      "Taxa are not removed across iterations to maintain consistent dimensions. \nDownstream analyses should handle zero-abundance taxa appropriately."
    )
  }

  # Data Sparsity
  cli::cli_h3("Data Sparsity")

  if (is_array) {
    iter_zeros <- apply(processed_data, 3, function(mat) sum(mat == 0))
    iter_total <- dim(processed_data)[1] * dim(processed_data)[2]

    cli::cli_ul()
    cli::cli_li("Rarefied matrix (across {.val {num_iter}} iterations):")
    cli::cli_ul(c(
      "Min: {.val {min(iter_zeros)}} zeros ({.val {round(min(iter_zeros) / iter_total * 100, 2)}}% sparsity) out of {.val {iter_total}} entries",
      "Max: {.val {max(iter_zeros)}} zeros ({.val {round(max(iter_zeros) / iter_total * 100, 2)}}% sparsity) out of {.val {iter_total}} entries",
      "Avg: {.val {round(mean(iter_zeros), 1)}} zeros ({.val {round(mean(iter_zeros) / iter_total * 100, 2)}}% sparsity) out of {.val {iter_total}} entries"
    ))
    cli::cli_end()
    cli::cli_end()
  } else {
    cli::cli_alert_info(
      "Returning list of data frames for each iteration."
    )

    rarefied_sparsity <- .sparsity_count(processed_data)
    # cli::cli_alert_info(
    #   "Rarefied matrix (across {.val {num_iter}} iterations):"
    # )
    cli::cli_ul()
    cli::cli_li("Rarefied matrix (across {.val {num_iter}} iterations):")
    cli::cli_ul(c(
      "Min: {.val {rarefied_sparsity$min_zeros}} zeros ({.val {rarefied_sparsity$min_sparsity}}% sparsity) out of {.val {rarefied_sparsity$min_total}} entries",
      "Max: {.val {rarefied_sparsity$max_zeros}} zeros ({.val {rarefied_sparsity$max_sparsity}}% sparsity) out of {.val {rarefied_sparsity$max_total}} entries",
      "Avg: {.val {round(rarefied_sparsity$avg_zeros, 1)}} zeros ({.val {rarefied_sparsity$avg_sparsity}}% sparsity) out of {.val {round(rarefied_sparsity$avg_total, 1)}} entries"
    ))
    cli::cli_end()
    cli::cli_end()
  }

  # Final Data Dimensions
  cli::cli_h3("Final Data Dimensions")

  if (is_array) {
    iter_nonzero_taxa <- apply(processed_data, 3, function(mat) {
      sum(colSums(mat) > 0)
    })
    n_samp <- dim(processed_data)[1]
    n_taxa <- dim(processed_data)[2]

    cli::cli_alert_success(
      "Output: {.val {num_iter}} iteration{?s} with {.val {n_samp}} sample{?s}"
    )
    cli::cli_ul()
    cli::cli_li("Array dimensions:")
    cli::cli_ul(c(
      "Samples: {.val {n_samp}}",
      "Taxa: {.val {n_taxa}}",
      "Iterations: {.val {num_iter}}"
    ))
    cli::cli_end()

    cli::cli_li("Non-zero taxa per iteration:")
    cli::cli_ul(c(
      "Min: {.val {min(iter_nonzero_taxa)}}",
      "Max: {.val {max(iter_nonzero_taxa)}}",
      "Avg: {.val {round(mean(iter_nonzero_taxa), 1)}}"
    ))
    cli::cli_end()
    cli::cli_end()
  } else {
    processed_data <- lapply(processed_data, function(df) {
      df[, colSums(df) > 0, drop = FALSE]
    })
    summary_stats <- list(
      min_taxa = min(vapply(processed_data, ncol, integer(1))),
      max_taxa = max(vapply(processed_data, ncol, integer(1))),
      avg_taxa = mean(vapply(processed_data, ncol, double(1))),
      min_samples = min(vapply(processed_data, nrow, integer(1))),
      max_samples = max(vapply(processed_data, nrow, integer(1)))
    )

    cli::cli_alert_success(
      "Output: {.val {num_iter}} iteration{?s} with {.val {length(unique_samples)}} unique sample{?s}"
    )
    cli::cli_ul()
    cli::cli_li("Samples per iteration:")
    cli::cli_ul(c(
      "Min: {.val {summary_stats$min_samples}}",
      "Max: {.val {summary_stats$max_samples}}"
    ))
    cli::cli_end()

    cli::cli_li("Non-zero taxa per iteration:")
    cli::cli_ul(c(
      "Min: {.val {summary_stats$min_taxa}}",
      "Max: {.val {summary_stats$max_taxa}}",
      "Avg: {.val {round(summary_stats$avg_taxa, 1)}}"
    ))
    cli::cli_end()
    cli::cli_end()
  }
}
