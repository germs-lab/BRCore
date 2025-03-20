#' Calculate Bray-Curtis Dissimilarity Between Sample Pairs
#'
#' This function calculates the Bray-Curtis dissimilarity between all pairs of samples in a matrix.
#' The dissimilarity is normalized by the total number of reads (`nReads`) to account for differences
#' in sequencing depth.
#'
#' @param matrix A numeric matrix or data frame where rows represent taxa (e.g., ASVs) and columns
#'   represent samples. The column names of the matrix are used to generate pairwise sample names.
#' @param nReads A numeric value representing the total number of reads used for normalization.
#'   Typically, this is the minimum or average sequencing depth across samples.
#' @return A list containing:
#'   \itemize{
#'     \item \code{values}: A numeric vector of Bray-Curtis dissimilarity values for each sample pair.
#'     \item \code{names}: A character vector of pairwise sample names (e.g., "Sample1-Sample2").
#'   }
#' @import utils
#' @export

calculate_bc <- function(matrix, nReads) {
  if (nrow(matrix) == 0) {
    cli::cli_alert_warning("{.arg matrix} is empty. Enter a non-empty matrix.")
    return(list(values = numeric(0), names = character(0)))
  }

  bc_values <- apply(combn(ncol(matrix), 2), 2, function(cols) {
    sum(abs(matrix[, cols[1]] - matrix[, cols[2]])) / (2 * nReads)
  })

  x_names <- apply(combn(ncol(matrix), 2), 2, function(cols) {
    paste(colnames(matrix)[cols], collapse = "-")
  })

  list(values = bc_values, names = x_names)
}
