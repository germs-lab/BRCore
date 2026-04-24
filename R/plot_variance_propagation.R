#' Variance propagation diagnostic for rarefaction
#'
#' This function evaluate the variance generated during multiple rarefaction by plotting comparing raw vs. rarefied alpha diversity metrics calculated at each iterations. It is possible to plot observed richness (q=0), Shannon diversity (q=1), or Simpson diversity (q=2) by setting the `q` parameter to "richness" or `q` = 0, "shannon" or `q` = 1, or "shannon" or `q` = 2. The plot is faceted by method (raw vs rarefied) and colored by a specified grouping variable from the sample data.
#'
#' @param physeq_obj Raw phyloseq object
#' @param rarefied Output from multi_rarefy(). Either a list of dataframes or and array.
#' @param q Hill number order (q = 0 for richness, q = 1 for Shannon, q = 2 for Simpson)
#' @param group_var A grouping variable to use gor grouping as in the sample_data()
#' @param group_color A color variable to use present in the sample_data()
#' @param convert_to_factor Logical. If \code{TRUE}, both \code{group_var} and
#'   \code{group_color} are coerced to \code{factor} before plotting, which is
#'   useful when those columns are numeric/continuous (e.g. dates, counts) but
#'   should be treated as discrete groups. When \code{TRUE} a discrete colour
#'   scale (\code{scale_color_viridis_d}) is used; otherwise the continuous
#'   scale (\code{scale_color_viridis_c}) is used. Default \code{FALSE}.
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
#' bcse_filt <- bcse |>
#' subset_samples(Crop %in% c("Poplar", "Switchgrass"))
#' bcse_rarefied_otutable_filt <-
#'  multi_rarefy(
#'        physeq_obj = bcse_filt,
#'        depth_level = 1000,
#'        num_iter = 100,
#'        .as = "list",
#'        set_seed = 7643
#'    )
#'
#'plot_variance_propagation(
#'    physeq_obj   = bcse_filt,
#'    rarefied = bcse_rarefied_otutable_filt,
#'    q        = 1,
#'    group_var = "Crop",
#'    group_color = "Plot"
#')
#' }
#'
#' @export
plot_variance_propagation <- function(
  physeq_obj,
  rarefied,
  q = 0,
  group_var,
  group_color,
  convert_to_factor = FALSE
) {
  # Validate inputs
  .phyloseq_class_check(physeq_obj)
  .validate_group_vars(physeq_obj, group_var, group_color)

  cli::cli_h1("Rarefaction Variance Propagation Visualization")

  # Extract raw OTU

  otu_raw <- .extract_otu_matrix(physeq_obj, samples_as_rows = TRUE)

  metadata <- data.frame(phyloseq::sample_data(physeq_obj))
  metadata$sample_id <- rownames(metadata)

  rare_samples <- .get_sample_ids(rarefied)

  common_samples <- intersect(rownames(otu_raw), rare_samples)

  otu_raw <- otu_raw[common_samples, , drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  .validate_group_vars(physeq_obj, group_var, group_color)

  cli::cli_alert_info("Hill number order selected, q= {.val {q}}")

  raw_vals <- .hill_number(otu_raw, q = q)

  raw_df <- data.frame(
    sample_id = names(raw_vals),
    value = raw_vals,
    method = "Raw"
  ) |>
    left_join(
      metadata[, c("sample_id", group_var, group_color), drop = FALSE],
      by = "sample_id"
    )

  n_iter <- .get_n_iterations(rarefied)

  cli::cli_alert_info(
    "Number of rarefaction iterations, n_iter= {.val {n_iter}}"
  )

  # replicate raw values (metadata already joined)
  raw_df <- raw_df[rep(seq_len(nrow(raw_df)), n_iter), ]

  iter_list <- .get_iter_list(rarefied)

  # Check rarefied depth
  rarefied_depths <- vapply(
    iter_list,
    function(mat) {
      rowSums(as.matrix(mat))[1]
    },
    double(1)
  )

  rare_df <- lapply(seq_along(iter_list), function(i) {
    mat <- as.matrix(iter_list[[i]])
    vals <- .hill_number(mat, q = q)
    data.frame(
      sample_id = names(vals),
      value = vals,
      method = "Rarefied"
    )
  }) |>
    bind_rows() |>
    left_join(
      metadata[, c("sample_id", group_var, group_color), drop = FALSE],
      by = "sample_id"
    )

  # Combine
  all_df <- bind_rows(
    mutate(raw_df, iter = NA_integer_),
    mutate(rare_df, iter = rep(seq_len(n_iter), each = length(common_samples)))
  ) |>
    mutate(method = factor(.data$method, levels = c("Raw", "Rarefied")))

  # Optionally coerce grouping columns to factor
  if (convert_to_factor) {
    all_df <- all_df |>
      mutate(across(all_of(c(group_var, group_color)), as.factor))
    cli::cli_alert_info(
      "Converted {.val {group_var}} and {.val {group_color}} to factor."
    )
  }

  # Choose color scale based on whether group_color is discrete
  color_scale <- if (
    convert_to_factor ||
      is.factor(all_df[[group_color]]) ||
      is.character(all_df[[group_color]])
  ) {
    ggplot2::scale_color_viridis_d(option = "turbo")
  } else {
    ggplot2::scale_color_viridis_c(option = "turbo")
  }

  p <- ggplot(
    all_df,
    aes(x = .data[[group_var]], y = .data$value, color = .data[[group_color]])
  ) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 0.8) +
    facet_wrap(~ .data$method) +
    .brcore_theme(axis_text_angle = 90) +
    color_scale +
    guides(
      color = guide_legend(override.aes = list(size = 5, shape = 15))
    ) +
    labs(
      x = group_var,
      title = "Raw vs Rarefied Alpha Diversity",
      subtitle = paste0(
        "Rarefied to depth of ",
        rarefied_depths[1],
        ", over ",
        n_iter,
        " iterations"
      ),
      y = if (is.numeric(q)) {
        switch(
          as.character(q),
          "0" = "Hill number (q = 0)",
          "1" = "Hill number (q = 1)",
          "2" = "Hill number (q = 2)",
          paste0("General Hill number (q = ", q, ")")
        )
      } else {
        switch(
          q,
          "richness" = "Observed Richness",
          "shannon" = "Shannon Diversity Index",
          "simpson" = "Simpson Diversity Index",
          paste0("Hill number ", q)
        )
      }
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
    q_val <- match.arg(q, INDICES)
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

  # Compute diversity based on q
  if (q_val == 0 || q_val == "richness") {
    # Richness: count non-zero species (same regardless of hill_form)
    if (length(dim(x)) > 1) {
      H <- rowSums(!is.na(x)) # Count non-NA values
    } else {
      H <- sum(!is.na(x))
    }
  } else if (q_val == 1 || q_val == "shannon") {
    # Shannon
    x_log <- x * log(x)
    if (length(dim(x)) > 1) {
      H <- -rowSums(x_log, na.rm = TRUE)
    } else {
      H <- -sum(x_log, na.rm = TRUE)
    }

    # Convert to Hill number if requested
    if (q_val == 1) {
      H <- exp(H)
    }
  } else if (q_val == 2 || q_val == "simpson") {
    # Simpson
    x_sq <- x^2
    if (length(dim(x)) > 1) {
      H <- rowSums(x_sq, na.rm = TRUE)
    } else {
      H <- sum(x_sq, na.rm = TRUE)
    }

    if (q_val == 2) {
      #Inverse Simpson concentration: 1 / sum(p^2)
      H <- 1 / H
    } else {
      # Simpson index: 1 - sum(p^2)
      H <- 1 - H
    }
  } else {
    # General Hill number not 0,1,2: (sum(p^q))^(1/(1-q))
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
#' @param physeq_obj A phyloseq object
#' @param group_var Character, grouping variable name
#' @param group_color Character, color variable name (optional)
#' @return Invisible NULL (errors if validation fails)
#' @keywords internal
#' @noRd
.validate_group_vars <- function(physeq_obj, group_var, group_color = NULL) {
  sample_df <- as(phyloseq::sample_data(physeq_obj), "data.frame")

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
