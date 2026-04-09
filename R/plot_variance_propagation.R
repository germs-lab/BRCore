#' Variance propagation diagnostic for rarefaction
#'
#' This function evaluate the variance generated during multiple rarefaction iterations.
#' It compares raw vs rarefied diversity metrics calculated at each iterations.
#'
#' @param physeq_obj Raw phyloseq object
#' @param rarefied Output from multi_rarefy(). Either a list of dataframes or and array.
#' @param q Hill number order (q = 0 for richness, q = 1 for Shannon, q = 2 for Simpson)
#' @param distance "bray" or "jaccard" (beta only)
#' @param group_var A grouping variable to use gor grouping as in the sample_data()
#' @param group_color A color variable to use present in the sample_data()
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table sample_data taxa_are_rows
#' @importFrom cli cli_h1 cli_h2 cli_alert_info cli_alert_warning cli_alert_success cli_alert_danger
#' @importFrom utils head
#' @importFrom vegan rrarefy
#' @importFrom dplyr left_join
#'
#' @return ggplot object comparing raw vs rarefied diversity distributions across iterations.
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' # Example using the bcse dataset, comparing hill q=1 between Poplar and Switchgrass plots
#' bcse_filt <- bcse %>%
#' subset_samples(Crop %in% c("Poplar", "Switchgrass"))
#' bcse_rarefied_otutable_filt <-
#'  multi_rarefy(
#'        physeq = bcse_filt,
#'        depth_level = 1000,
#'        num_iter = 100,
#'        .as_array = FALSE,
#'        set_seed = 7643
#'    )
#'
#'plot_variance_propagation(
#'    physeq   = bcse_filt,
#'    rarefied = bcse_rarefied_otutable_filt,
#'    q        = 1,
#'    group_var = "Crop",
#'    group_color = "Plot"
#') + scale_color_viridis_d(option = "turbo")
#' }
#'
#' @export
plot_variance_propagation <- function(
  physeq_obj,
  rarefied,
  q = 0,
  group_var,
  group_color
) {
  cli::cli_h1("Rarefaction Variance Propagation Visualization")

  # Extract raw OTU

  otu_raw <- .extract_otu_matrix(physeq_obj, samples_as_rows = TRUE)

  metadata <- data.frame(phyloseq::sample_data(physeq_obj))
  metadata$sample_id <- rownames(metadata)

  rare_samples <- .get_sample_ids(rarefied)

  common_samples <- intersect(rownames(otu_raw), rare_samples)

  otu_raw <- otu_raw[common_samples, , drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  cli::cli_alert_info("Hill number order selected, q= {.val {q}}")

  raw_vals <- .hill_number(otu_raw, q = q)

  raw_df <- data.frame(
    sample_id = names(raw_vals),
    value = raw_vals,
    method = "Raw"
  )

  n_iter <- .get_n_iterations(rarefied)

  cli::cli_alert_info(
    "Number of rarefaction iterations, n_iter= {.val {n_iter}}"
  )

  # replicate raw values
  raw_df <- raw_df[rep(seq_len(nrow(raw_df)), n_iter), ]

  iter_list <- .get_iter_list(rarefied)

  rare_df <- lapply(seq_along(iter_list), function(i) {
    mat <- as.matrix(iter_list[[i]])

    vals <- .hill_number(mat, q = q)

    data.frame(
      sample_id = names(vals),
      value = vals,
      method = "Rarefied"
    )
  }) %>%
    bind_rows()

  # Combine
  all_df <- bind_rows(
    mutate(raw_df, iter = NA_integer_),
    mutate(rare_df, iter = rep(seq_len(n_iter), each = length(common_samples)))
  ) %>%
    left_join(
      metadata[, c("sample_id", group_var, group_color)],
      by = "sample_id"
    )

  p <- ggplot(
    all_df,
    aes(x = .data[[group_var]], y = value, color = .data[[group_color]])
  ) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 0.8) +
    facet_wrap(~method) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.background = element_blank(),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 5, shape = 15)
      )
    ) +
    labs(
      x = group_var,
      y = paste0("Hill number (q = ", q, ")"),
      title = "Raw vs Rarefied Alpha Diversity"
    )

  cli::cli_alert_info("Comparison plot generated!")
  return(p)
}


#' Compute Hill numbers for samples
#'
#' Based on vegan::diversity() function, but modified to
#' handle both vector and matrix input, and to support named
#' indices for common q values.
#' @param x Numeric vector or matrix of species abundances.
#'   If matrix, samples are rows and taxa are columns.
#' @param q Hill number order. Can be:
#'   - "richness" (q=0): Species richness
#'   - "shannon" (q=1): Exponential Shannon entropy
#'   - "simpson" (q=2): Inverse Simpson concentration
#'   - Numeric value: General Hill number of order q
#' @return Numeric vector of Hill numbers (one per sample)
#' @keywords internal
#' @noRd
.hill_number <- function(x, q = "richness") {
  # Named indices for common cases
  INDICES <- c("richness", "shannon", "simpson")

  # Handle numeric q or match named index
  if (is.numeric(q)) {
    q_val <- q
  } else {
    index <- match.arg(q, INDICES)
    q_val <- switch(index, richness = 0, shannon = 1, simpson = 2)
  }

  # Handle matrix (rows = samples) or vector input
  if (length(dim(x)) > 1) {
    # Matrix: rows are samples, columns are taxa
    total <- rowSums(x)
    x <- x / total # Proportions (vectorized)
    x[x == 0] <- NA # Set zeros to NA (ignored in rowSums)
  } else {
    # Vector: single sample
    total <- sum(x)
    x <- x / total
    x[x == 0] <- NA
  }

  # Compute Hill number based on q
  if (q_val == 0) {
    # Richness: count non-zero species
    if (length(dim(x)) > 1) {
      H <- rowSums(!is.na(x))
    } else {
      H <- sum(!is.na(x))
    }
  } else if (q_val == 1) {
    # Shannon: exp(-sum(p * log(p)))
    x_log <- x * log(x)
    if (length(dim(x)) > 1) {
      H <- exp(-rowSums(x_log, na.rm = TRUE))
    } else {
      H <- exp(-sum(x_log, na.rm = TRUE))
    }
  } else {
    # General Hill number: (sum(p^q))^(1/(1-q))
    x_pow <- x^q_val
    if (length(dim(x)) > 1) {
      H <- rowSums(x_pow, na.rm = TRUE)^(1 / (1 - q_val))
    } else {
      H <- sum(x_pow, na.rm = TRUE)^(1 / (1 - q_val))
    }
  }

  # Handle NA totals (samples with zero abundance)
  if (any(NAS <- is.na(total))) {
    H[NAS] <- NA
  }

  H
}

#' Validate that group variables exist in sample data
#'
#' @param physeq A phyloseq object
#' @param group_var Character, grouping variable name
#' @param group_color Character, color variable name (optional)
#' @return Invisible NULL (errors if validation fails)
#' @keywords internal
#' @noRd
.validate_group_vars <- function(physeq, group_var, group_color = NULL) {
  sample_df <- as(phyloseq::sample_data(physeq), "data.frame")

  if (!group_var %in% colnames(sample_df)) {
    cli::cli_abort(
      "{.arg group_var} '{group_var}' not found in sample_data()."
    )
  }

  if (!is.null(group_color) && !group_color %in% colnames(sample_df)) {
    cli::cli_abort(
      "{.arg group_color} '{group_color}' not found in sample_data()."
    )
  }

  invisible(NULL)
}
