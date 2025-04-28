#' Calculate NMDS for ASV Communities
#'
#' This function performs Hellinger transformation, calculates Bray-Curtis distances,
#' and computes NMDS ordination for ASV communities. It includes error handling and
#' checks for NA/NaN values.
#'
#' @param asv_matrix A matrix or data frame of ASV counts (samples x ASVs).
#' @param physeq A phyloseq object containing sample metadata.
#' @param ncores Number of cores for parallel processing. Default is `parallel::detectCores()`.
#' @param k Number of dimensions for NMDS. Default is 2.
#' @param maxit Maximum number of iterations for NMDS. Default is 999.
#' @param trymax Maximum number of random starts for NMDS. Default is 100.
#' @param previous.best A previous NMDS result to use as a starting point.
#' @return A list containing:
#'   \itemize{
#'     \item \code{nmds_scores}: NMDS scores for samples.
#'     \item \code{nmds_df}: NMDS scores merged with sample metadata.
#'     \item \code{ordi}: NMDS ordination object.
#'   }
#' @import phyloseq
#' @import dplyr
#' @import vegan
#' @export
brc_nmds <- function(asv_matrix,
                           physeq,
                           ncores = parallel::detectCores(),
                           k = 2,
                           maxit = 999,
                           trymax = 100,
                           previous.best = NULL) {
    set.seed(54641)
    
    # Validate inputs
    if (!is.matrix(asv_matrix) && !is.data.frame(asv_matrix)) {
        cli::cli_abort("{.arg asv_matrix} must be a matrix or data frame.")
    }
    
    if (!inherits(physeq, "phyloseq")) {
        cli::cli_abort("{.arg physeq} must be a phyloseq object.")
    }
    
    if (any(is.na(asv_matrix)) || any(is.nan(asv_matrix))) {
        cli::cli_abort("{.arg asv_matrix} contains NA or NaN values.")
    }
    
    # Hellinger transformation
    hell_matrix <- decostand(t(asv_matrix), MARGIN = 1, method = "hellinger")
    
    # Remove columns that sum to 0
    hell_matrix <- hell_matrix %>%
        as.data.frame() %>%
        select(where( ~ is.numeric(.) && sum(.) > 0)) %>%
        as.matrix()
    
    # Check for empty matrix after filtering
    if (ncol(hell_matrix) == 0) {
        cli::cli_abort("No valid columns remaining after filtering.")
    }
    
    # Calculate Bray-Curtis distances
    asv_dist <- vegdist(
        t(hell_matrix),
        method = "bray",
        upper = FALSE,
        binary = FALSE,
        na.rm = TRUE
    )
    
    
    # Perform NMDS
    ordi <- metaMDS(
        as.matrix(asv_dist),
        distance = "bray",
        display = "sites",
        noshare = TRUE,
        autotransform = FALSE,
        wascores = TRUE,
        tidy = TRUE,
        k = k,
        maxit = 999,
        trymax = trymax,
        parallel = ncores,
        previous.best = previous.best
    )
    
    
    # Check NMDS convergence
    if (ordi$converged == FALSE) {
        cli::cli_alert_warning("NMDS did not converge. Consider increasing {.arg trymax}.")
    }
    
    # Add species scores
    vegan::sppscores(ordi) <- t(hell_matrix)
    
    # Extract NMDS scores
    nmds_scores <- as.data.frame(vegan::scores(ordi)$sites)
    nmds_scores <- nmds_scores %>%
        rownames_to_column(var = "unique_id")
    
    # Merge with sample metadata
    sample_df <- physeq@sam_data %>%
        data.frame() %>%
        rownames_to_column(var = "unique_id")
    
    nmds_df <- right_join(sample_df, nmds_scores, by = "unique_id")
    
    # Return results
    return(list(
        nmds_scores = nmds_scores,
        nmds_df = nmds_df,
        ordi = ordi
    ))
}
