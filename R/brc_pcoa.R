#' Calculate Principal Coordinate Analysis (PCoA) using the weighted classical multidimensional scaling (wcmdscale) method.

#' @param asv_matrix A matrix, data frame, or distance object containing the ASV (Amplicon Sequence Variant) data 
#'   to be used for Principal Coordinate Analysis.
#' @param physeq A phyloseq object containing sample metadata that will be incorporated into the PCoA results.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing. 
#'   Default is all available cores as detected by parallel::detectCores().
#' @param k Integer specifying the number of dimensions (axes) to calculate in the PCoA. Default is 2.
#' @param eig Logical indicating whether eigenvalues should be returned. Default is TRUE.
#' @param add Logical indicating whether to add a constant to the non-diagonal dissimilarities such that all 
#'   eigenvalues are non-negative. Default is FALSE.
#' @param x.ret Logical indicating whether the rotated configuration should be returned. Default is FALSE.
#' 
#' #' @export

brc_pcoa <- function(
  asv_matrix,
  physeq,
  ncores = parallel::detectCores(),
  k = 2,
  eig = TRUE,
  add = FALSE,
  x.ret = FALSE
) {
  set.seed(54641)

  # Validate inputs
  if (
    !is.matrix(asv_matrix) &&
      !is.data.frame(asv_matrix) &&
      !"dist" %in% class(asv_matrix)
  ) {
    cli::cli_abort(
      "{.arg asv_matrix} must be a matrix, distance or data frame."
    )
  }

  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort("{.arg physeq} must be a phyloseq object.")
  }

  if (any(is.na(asv_matrix)) || any(is.nan(asv_matrix))) {
    cli::cli_abort("{.arg asv_matrix} contains NA or NaN values.")
  }

  # Perform PCoA
  ordi <- wcmdscale(asv_matrix, k = k, eig = eig, add = add, x.ret = x.ret)

  # Extract PCoA scores and add metadata
  pcoa_df <- data.frame(
    Dim1 = ordi$points[, 1], # First PCoA axis
    Dim2 = ordi$points[, 2], # Second PCoA axis
    brc = factor(physeq@sam_data$brc), # Convert 'brc' to a factor
    crop = factor(physeq@sam_data$crop) # Convert 'crop' to a factor
  )

  # Return results
  return(list(
    pcoa_df = pcoa_df,
    ordi = ordi,
    distance_matrix = asv_matrix # Added distance matrix to output
  ))
}
