# R/utils-otu-extraction.R

#' Check if input is a phyloseq object
#' @param physeq_obj An object to check
#' @return Invisible NULL (errors if not a
#' phyloseq object)
#' @keywords internal
#' @noRd

.phyloseq_class_check <- function(physeq_obj) {
  if (inherits(physeq_obj, "phyloseq")) {
    cli::cli_alert_success("Input phyloseq object is valid!")
  }
  if (!inherits(physeq_obj, "phyloseq")) {
    cli::cli_abort(
      "Input must be a 'phyloseq' object, not a {.cls {class(physeq_obj)}}"
    )
  }
}


#' Extract OTU matrix from phyloseq object
#'
#' @param physeq_obj A phyloseq object
#' @return A matrix with samples as rows and taxa as columns
#' @keywords internal
#' @noRd
.extract_otu_matrix <- function(physeq_obj, samples_as_rows = TRUE) {
  suppressMessages(.phyloseq_class_check(physeq_obj))

  otu_mat <- as(phyloseq::otu_table(physeq_obj), "matrix")

  if (phyloseq::taxa_are_rows(physeq_obj) && samples_as_rows) {
    otu_mat <- t(otu_mat)
  } else if (!phyloseq::taxa_are_rows(physeq_obj) && !samples_as_rows) {
    otu_mat <- t(otu_mat)
  }

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
  } else if (is.list(rarefied)) {
    return(length(rarefied))
  } else {
    cli::cli_abort(
      "{.arg rarefied} must be a list or 3D array from multi_rarefy()."
    )
  }
}

#' Convert rarefied object to a list of matrices
#'
#' @param x A 3D array or list of matrices from `multi_rarefy()`
#' @return A list of matrices, one per iteration
#' @keywords internal
#' @noRd
.get_iter_list <- function(x) {
  if (is.array(x) && length(dim(x)) == 3) {
    lapply(seq_len(dim(x)[3]), function(i) x[,, i])
  } else if (is.list(x)) {
    x
  } else {
    list(as.matrix(x))
  }
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
  } else if (is.list(rarefied)) {
    sample_ids <- rownames(rarefied[[1]])

    return(sample_ids)
  } else {
    cli::cli_abort(
      "{.arg rarefied} must be a list or 3D array from multi_rarefy()."
    )
  }
}
