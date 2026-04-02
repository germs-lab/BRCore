#' Run multiple rarefaction for microbiome count tables
#'
#' This function performs multiple rarefaction on a `phyloseq` object by randomly
#' sub-sampling OTUs/ASVs within samples without replacement. The process is
#' repeated for a specified number of iterations, and the results are averaged.
#' Samples with fewer OTUs/ASVs than the specified `depth_level` are discarded.
#'
#' @param physeq A `phyloseq` object containing an OTU/ASV table.
#' @param depth_level An integer specifying the sequencing depth (number of
#'   OTUs/ASVs) to which samples should be rarefied.
#' @param num_iter An integer specifying the number of iterations to perform
#'   for rarefaction.
#' @param .summarize A logical indicating whether to summarize the results by
#'   averaging across iterations (default = TRUE). If FALSE, the function will
#'   return a named list of data frames (e.g. iter_1), one for each iteration, without
#' averaging.
#' @param threads Number of threads (default = `get_available_cores()`).
#' @param set_seed An optional integer to set the random seed for reproducibility (default = NULL).
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
#'   physeq = bcse,
#'   depth_level = 1000,
#'   num_iter = 100,
#'   .summarize = TRUE,
#'   threads = 2,
#'   set_seed = 7642
#' )
#'
#' rowSums(otu_table_rare)
#' }
#'
#' @export
multi_rarefy <- function(
  physeq,
  depth_level,
  num_iter = 100,
  .summarize = TRUE,
  threads = get_available_cores(),
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

  if (!inherits(physeq, "phyloseq")) {
    cli::cli_alert_danger(
      "Input must be a phyloseq object, not a {.cls {class(physeq)}}"
    )
    stop("Input must be a phyloseq object")
  }

  if (!is.numeric(depth_level) || depth_level <= 0) {
    cli::cli_alert_danger("depth_level must be a positive number")
    stop("Invalid depth_level")
  }

  if (!is.numeric(num_iter) || num_iter <= 0) {
    cli::cli_alert_danger("num_iter must be a positive integer")
    stop("Invalid num_iter")
  }

  if (!is.numeric(threads) || threads <= 0) {
    cli::cli_alert_danger("threads must be a positive integer")
    stop("Invalid threads")
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
  otu_mat <- as.matrix(otu_table(physeq, ))

  if (phyloseq::taxa_are_rows(physeq)) {
    otu_mat <- t(otu_mat) # We want samples as rows
  }

  dataframe <- as.data.frame(otu_mat)

  # Input parameter checks ----
  cli::cli_alert_info(
    "Input (matrix/df dim): {.val {nrow(dataframe)}} samples x {.val {ncol(dataframe)}} taxa"
  )
  cli::cli_alert_info("Rarefaction depth: {.val {depth_level}}")
  cli::cli_alert_info("Iterations: {.val {num_iter}}")

  cli::cli_alert_info("taxa_are_rows: {phyloseq::taxa_are_rows(physeq)}")

  cli::cli_alert_info(
    "OTU matrix/df rownames head: {paste(head(rownames(dataframe)), collapse = ', ')}"
  )
  cli::cli_alert_info(
    "OTU matrix/df colnames head: {paste(head(colnames(dataframe)), collapse = ', ')}"
  )

  cli::cli_alert_info(
    "Row sums summary: Min={min(rowSums(dataframe))}, Max={max(rowSums(dataframe))}, Median={median(rowSums(dataframe))}"
  )

  # Parallel setup ----
  threads <- min(threads, availableCores())
  cl <- makeCluster(threads)
  on.exit(stopCluster(cl), add = TRUE)

  # Export needed objects/packages to workers
  clusterExport(
    cl,
    varlist = c(
      "dataframe",
      "depth_level",
      "iteration_seeds"
    ),
    envir = environment()
  )
  clusterEvalQ(cl, library(vegan))

  # Run rarefactions in parallel ----
  cli::cli_alert_info("Running rarefaction...")

  com_iter_list <- vector(mode = "list", length = num_iter)
  com_iter_list <- parLapply(cl, 1:num_iter, function(i) {
    if (!is.na(iteration_seeds[i])) {
      set.seed(iteration_seeds[i])
    }

    df <- vegan::rrarefy(
      dataframe,
      sample = depth_level
    )

    # Summarization metadata
    df <- as.data.frame(df) |>
      rownames_to_column("sample_id") |>
      mutate(iter = i)
  })

  # Aggregate results and removal of zero-abundance taxa ----
  n_samples_before <- nrow(dataframe)

  if (.summarize) {
    com_iter_df <- do.call(rbind, com_iter_list) |>
      as.data.frame()

    processed_data <- com_iter_df |>
      group_by(sample_id) |>
      summarise(across(everything(), mean), .groups = "drop") |>
      select(!iter) |>
      filter(
        dplyr::near(
          # We use near() to account for floating-point precision issues that can arise when averaging counts
          rowSums(across(where(is.numeric))),
          depth_level,
          tol = 1e-6 # Accept values within ±0.5 of depth_level
        )
      ) |>
      column_to_rownames("sample_id")

    # Cleaning
    processed_data <- processed_data[, colSums(processed_data) > 0]
  }

  if (!.summarize) {
    processed_data <- lapply(seq_along(com_iter_list), function(i) {
      df <- com_iter_list[[i]] |>
        mutate(
          sample_id = paste0(sample_id, "_iter_", iter)
        ) |>
        select(!c(iter)) |>
        filter(
          dplyr::near(
            rowSums(across(where(is.numeric))),
            depth_level,
            tol = 1e-6
          )
        ) |>
        column_to_rownames("sample_id")
    })
    names(processed_data) <- paste0("iter_", seq_len(num_iter))

    # Remove zero-abundance taxa for each iteration
    processed_data <- lapply(processed_data, function(df) df[, colSums(df) > 0])

    cli::cli_alert_info(
      "Returning list of data frames for each iteration (unsummarized)."
    )
  }

  if (.summarize) {
    n_taxa_before <- ncol(processed_data)
    n_taxa_after <- ncol(processed_data)

    n_samples_after <- nrow(processed_data)
    n_samples_removed <- n_samples_before - n_samples_after
    removed_samples <- setdiff(rownames(dataframe), rownames(processed_data))
  } else {
    n_taxa_before <- ncol(processed_data[[1]])
    n_taxa_after <- min(sapply(processed_data, ncol))

    n_samples_after <- nrow(processed_data[[1]])
    n_samples_removed <- n_samples_before - n_samples_after
    removed_samples <- setdiff(
      rownames(dataframe),
      rownames(processed_data[[1]])
    )

    unique_samples <- length(unique(sub(
      "_iter_.*",
      "",
      rownames(processed_data[[1]])
    )))
    n_samples_removed <- n_samples_before - unique_samples
    original_sample_ids <- sub("_iter_.*", "", rownames(processed_data[[1]])) |>
      unique()
    removed_samples <- setdiff(rownames(dataframe), original_sample_ids)
  }

  n_taxa_removed <- n_taxa_before - n_taxa_after

  .report_rarefaction_results(
    n_samples_removed = n_samples_removed,
    n_taxa_removed = n_taxa_removed,
    removed_samples = removed_samples,
    dataframe = dataframe,
    processed_data = processed_data,
    depth_level = depth_level,
    .summarize = .summarize,
    unique_samples = if (.summarize) NULL else unique_samples,
    num_iter = if (.summarize) NULL else num_iter
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
    min_zeros <- min(sapply(dataframe, function(df) sum(df == 0)))
    min_total <- min(sapply(dataframe, function(df) nrow(df) * ncol(df)))

    max_zeros <- max(sapply(dataframe, function(df) sum(df == 0)))
    max_total <- max(sapply(dataframe, function(df) nrow(df) * ncol(df)))

    avg_zeros <- mean(sapply(dataframe, function(df) sum(df == 0)))
    avg_total <- mean(sapply(dataframe, function(df) nrow(df) * ncol(df)))
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
#' @param .summarize Whether the data was summarized (TRUE) or not (FALSE).
#' @param unique_samples The number of unique samples in the processed data (only relevant if .summarize = FALSE).
#' @param num_iter The number of iterations performed (only relevant if .summarize = FALSE).
#' @return None. This function is used for side-effect reporting via cli messages.
#' @noRd
#' @keywords internal
#'
.report_rarefaction_results <- function(
  n_samples_removed,
  n_taxa_removed,
  removed_samples,
  dataframe,
  processed_data,
  depth_level,
  .summarize = TRUE,
  unique_samples = NULL,
  num_iter = NULL
) {
  cli::cli_h2("Rarefaction Results")

  # Sample Removal
  cli::cli_h3("Sample Removal")
  if (n_samples_removed > 0) {
    cli::cli_alert_warning(
      "{.val {n_samples_removed}} sample{?s} removed (depth < {.val {depth_level}})"
    )
    cli::cli_alert_warning(
      "Samples removed: {.val {paste(removed_samples, collapse = ', ')}}"
    )
  }

  # Taxa Removal
  cli::cli_h3("Taxa Removal")
  if (n_taxa_removed > 0) {
    cli::cli_alert_warning(
      "{.val {n_taxa_removed}} taxa removed (zero abundance)"
    )
  }

  # Data Sparsity
  cli::cli_h3("Data Sparsity")
  original_sparsity <- .sparsity_count(dataframe)
  cli::cli_alert_info(
    "Original matrix: {.val {original_sparsity$zeros}} zeros ({.val {original_sparsity$sparsity}}% sparsity) out of {.val {original_sparsity$total}} entries"
  )

  rarefied_sparsity <- .sparsity_count(processed_data)
  if (!.summarize) {
    cli::cli_alert_info(
      "Rarefied matrix: (across {.val {num_iter}} iterations):"
    )
    cli::cli_ul(c(
      "Min: {.val {rarefied_sparsity$min_zeros}} zeros ({.val {rarefied_sparsity$min_sparsity}}% sparsity) out of {.val {rarefied_sparsity$min_total}} entries",
      "Max: {.val {rarefied_sparsity$max_zeros}} zeros ({.val {rarefied_sparsity$max_sparsity}}% sparsity) out of {.val {rarefied_sparsity$max_total}} entries",
      "Avg: {.val {round(rarefied_sparsity$avg_zeros, 1)}} zeros ({.val {rarefied_sparsity$avg_sparsity}}% sparsity) out of {.val {round(rarefied_sparsity$avg_total, 1)}} entries"
    ))
  }

  # Final Data Dimensions
  cli::cli_h3("Final Data Dimensions")
  if (.summarize) {
    cli::cli_alert_success(
      "Output: {.val {nrow(processed_data)}} samples x {.val {ncol(processed_data)}} taxa"
    )
  } else {
    summary_stats <- list(
      min_taxa = min(sapply(processed_data, ncol)),
      max_taxa = max(sapply(processed_data, ncol)),
      avg_taxa = mean(sapply(processed_data, ncol)),
      min_samples = min(sapply(processed_data, nrow)),
      max_samples = max(sapply(processed_data, nrow))
    )

    cli::cli_alert_success(
      "Output: {.val {num_iter}} iteration{?s} with {.val {unique_samples}} unique sample{?s}"
    )
    cli::cli_ul()
    cli::cli_li("Samples per iteration:")
    cli::cli_ul(c(
      "Min: {.val {summary_stats$min_samples}}",
      "Max: {.val {summary_stats$max_samples}}"
    ))
    cli::cli_end()

    cli::cli_li("Taxa per iteration:")
    cli::cli_ul(c(
      "Min: {.val {summary_stats$min_taxa}}",
      "Max: {.val {summary_stats$max_taxa}}",
      "Avg: {.val {round(summary_stats$avg_taxa, 1)}}"
    ))
    cli::cli_end()

    cli::cli_end()
  }
}
