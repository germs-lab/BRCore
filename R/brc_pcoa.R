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
