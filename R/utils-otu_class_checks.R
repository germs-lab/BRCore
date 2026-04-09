# R/utils-otu-extraction.R

#' Check if input is a phyloseq object
#' @param physeq An object to check
#' @return Invisible NULL (errors if not a
#' phyloseq object)
#' @keywords internal
#' @noRd

.phyloseq_class_check <- function(physeq) {
  if (inherits(physeq, "phyloseq")) {
    cli::cli_alert_success("Input phyloseq object is valid!")
  }
  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort(
      "Input must be a 'phyloseq' object, not a {.cls {class(physeq)}}"
    )
  }
}


#' Extract OTU matrix from phyloseq object
#'
#' @param physeq A phyloseq object
#' @return A matrix with samples as rows and taxa as columns
#' @keywords internal
#' @noRd
.extract_otu_matrix <- function(physeq, samples_as_rows = TRUE) {
  .phyloseq_class_check(physeq)

  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")

  if (phyloseq::taxa_are_rows(physeq) && samples_as_rows) {
    otu_mat <- t(otu_mat)
  } else if (!phyloseq::taxa_are_rows(physeq) && !samples_as_rows) {
    otu_mat <- t(otu_mat)
  }
  # if (phyloseq::taxa_are_rows(physeq)) {
  #   otu_mat <- t(otu_mat)
  # }

  otu_mat
}

#' Extract OTU matrix from phyloseq object or rarefied list/array
#'
#' @param rarefied A phyloseq object, list of data frames, or 3D array
#' @param iteration Which iteration to extract (for list/array input)
#' @return A matrix with samples as rows and taxa as columns
#' @keywords internal
#' @noRd
.extract_rarefied_matrix <- function(rarefied, iteration = NULL) {
  # Handle phyloseq object
  if (inherits(rarefied, "phyloseq")) {
    return(.extract_otu_matrix(rarefied))
  }

  # Handle 3D array
  if (is.array(rarefied) && length(dim(rarefied)) == 3) {
    if (is.null(iteration)) {
      cli::cli_abort("{.arg iteration} must be specified for array input.")
    }
    if (iteration < 1 || iteration > dim(rarefied)[3]) {
      cli::cli_abort(
        "{.arg iteration} must be between 1 and {dim(rarefied)[3]}."
      )
    }
    return(rarefied[,, iteration])
  }

  # Handle list of data frames/matrices
  if (is.list(rarefied)) {
    if (is.null(iteration)) {
      cli::cli_abort("{.arg iteration} must be specified for list input.")
    }
    if (iteration < 1 || iteration > length(rarefied)) {
      cli::cli_abort(
        "{.arg iteration} must be between 1 and {length(rarefied)}."
      )
    }
    return(as.matrix(rarefied[[iteration]]))
  }

  cli::cli_abort(
    "{.arg rarefied} must be a phyloseq object, list, or 3D array."
  )
}

#' Get number of iterations from rarefied object
#'
#' @param rarefied A list of data frames or 3D array from multi_rarefy()
#' @return Integer number of iterations
#' @keywords internal
#' @noRd
.get_n_iterations <- function(rarefied) {
  if (is.array(rarefied) && length(dim(rarefied)) == 3) {
    return(dim(rarefied)[3])
  }

  if (is.list(rarefied)) {
    return(length(rarefied))
  }

  cli::cli_abort(
    "{.arg rarefied} must be a list or 3D array from multi_rarefy()."
  )
}

#' Extract sample IDs from rarefied object
#'
#' @param rarefied A list of data frames or 3D array from multi_rarefy()
#' @param remove_suffix Logical, remove "_iter_N" suffix from list rownames
#' @return Character vector of sample IDs
#' @keywords internal
#' @noRd
.get_sample_ids <- function(rarefied, remove_suffix = TRUE) {
  if (is.array(rarefied) && length(dim(rarefied)) == 3) {
    return(dimnames(rarefied)[[1]])
  }

  if (is.list(rarefied)) {
    sample_ids <- rownames(rarefied[[1]])

    if (remove_suffix) {
      # Remove "_iter_N" suffix
      sample_ids <- sub("_iter_\\d+$", "", sample_ids)
      sample_ids <- unique(sample_ids)
    }

    return(sample_ids)
  }

  cli::cli_abort(
    "{.arg rarefied} must be a list or 3D array from multi_rarefy()."
  )
}

#' Compute Hill number for a single sample
#'
#' @param x Numeric vector of species abundances
#' @param q Hill number order (0 = richness, 1 = Shannon, 2 = Simpson)
#' @return Numeric Hill number value
#' @keywords internal
#' @noRd
.hill_number <- function(x, q = 0) {
  x <- x[x > 0]

  if (length(x) == 0) {
    return(0)
  }

  p <- x / sum(x)

  if (q == 0) {
    return(length(x))
  } else if (q == 1) {
    return(exp(-sum(p * log(p))))
  } else {
    return((sum(p^q))^(1 / (1 - q)))
  }
}

#' Compute Hill numbers for all samples in a matrix
#'
#' @param otu_mat Matrix with samples as rows, taxa as columns
#' @param q Hill number order (0 = richness, 1 = Shannon, 2 = Simpson)
#' @return Numeric vector of Hill numbers (one per sample)
#' @keywords internal
#' @noRd
.compute_hill_numbers <- function(otu_mat, q = 0) {
  apply(otu_mat, 1, .hill_number, q = q)
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
