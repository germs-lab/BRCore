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
#' @param threads Number of threads (default = 4).
#' @param set_seed An optional integer to set the random seed for reproducibility (default = NULL).
#'
#' @return A data frame with taxa as rows and samples as columns. The values
#'   represent the average sequence counts calculated across all iterations.
#'   Samples with less than `depth_level` sequences are discarded.
#'
#'
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom dplyr bind_rows group_by summarise across everything filter near
#' @importFrom dplyr where
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table
#' @importFrom cli cli_h1 cli_h2 cli_alert_info cli_alert_warning cli_alert_success cli_alert_danger
#' @importFrom utils head
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
  dataframe <- as.data.frame(
    as.matrix(t(otu_table(physeq, taxa_are_rows = TRUE)))
  )

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
      ".single_rarefy",
      "iteration_seeds"
    ),
    envir = environment()
  )

  # Run rarefactions in parallel ----
  cli::cli_alert_info("Running rarefaction...")

  com_iter <- parLapply(cl, 1:num_iter, function(i) {
    if (!is.na(iteration_seeds[i])) {
      set.seed(iteration_seeds[i])
    }
    df <- .single_rarefy(dataframe, sample_size = depth_level)
    df <- as.data.frame(df)
    rownames_to_column(df, "sample_id")
  })

  # Aggregate results ----
  n_samples_before <- nrow(dataframe)

  mean_data <- bind_rows(com_iter) |>
    group_by(sample_id) |>
    summarise(across(everything(), mean), .groups = "drop") |>
    filter(
      dplyr::near(
        # We use near() to account for floating-point precision issues that can arise when averaging counts
        rowSums(across(where(is.numeric))),
        depth_level,
        tol = 1e-6 # Accept values within ±0.5 of depth_level
      )
    ) |>
    column_to_rownames("sample_id")

  # Remove ASVs/OTUs with zero total abundance
  n_taxa_before <- ncol(mean_data)
  mean_data <- mean_data[, colSums(mean_data) > 0]
  n_taxa_after <- ncol(mean_data)

  # Report results ---
  n_samples_after <- nrow(mean_data)
  n_samples_removed <- n_samples_before - n_samples_after
  n_taxa_removed <- n_taxa_before - n_taxa_after

  cli::cli_h2("Rarefaction Results")
  if (n_samples_removed > 0) {
    cli::cli_alert_warning(
      "{.val {n_samples_removed}} sample{?s} removed (depth < {.val {depth_level}})"
    )
    removed_samples <- setdiff(rownames(dataframe), rownames(mean_data))
    cli::cli_alert_warning(
      "Samples removed: {.val {paste(removed_samples, collapse = ', ')}}"
    )
  }

  if (n_taxa_removed > 0) {
    cli::cli_alert_info(
      "{.val {n_taxa_removed}} taxa removed (zero abundance)"
    )
  }

  cli::cli_alert_success(
    "Output: {.val {nrow(mean_data)}} samples x {.val {ncol(mean_data)}} taxa"
  )

  return(mean_data)
}


#' Single rarefaction iteration
#'
#' @param x A matrix with samples as rows and taxa as columns.
#' @param sample_size The number of sequences to sample.
#'
#' @return A matrix with rarefied counts.
#'
#' @noRd
#' @keywords internal
.single_rarefy <- function(x, sample_size) {
  x <- as.matrix(x)

  out <- t(apply(x, 1, function(row) {
    total <- sum(row)
    if (total < sample_size) {
      return(rep(NA_integer_, length(row)))
    }

    sampled <- sample(
      rep(seq_along(row), times = row),
      size = sample_size,
      replace = FALSE
    )
    tabulate(sampled, nbins = length(row))
  }))

  colnames(out) <- colnames(x)
  rownames(out) <- rownames(x)
  out
}
