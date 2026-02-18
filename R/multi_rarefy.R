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
#' @details
#' To ensure reproducibility across different operating systems (macOS, Linux, Windows),
#' this function explicitly sets the RNG algorithm to "Mersenne-Twister" with "Inversion"
#' normal kind and "Rejection" sample kind before setting the seed. For parallel execution,
#' each iteration is assigned a pre-computed seed derived from the main seed, ensuring that
#' results are identical regardless of which worker processes which iteration. This approach
#' guarantees cross-platform reproducibility and thread-count independence.
#'
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom dplyr  bind_rows group_by summarise across everything filter
#' @importFrom dplyr where
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom phyloseq otu_table
#' @importFrom vegan rrarefy
#' @importFrom cli cli_text cli_warn
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' library(BRCore)
#' data("bcse", package = "BRCore")
#'
#' # Example rarefaction (single iteration, single core to keep examples fast)
#' otu_table_rare <- multi_rarefy(
#'    physeq = bcse,
#'    depth_level = 1000,
#'    num_iter = 100,
#'    threads = 2,
#'    set_seed = 7642
#')
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
    cli::cli_text("\nSeed used: {set_seed}\n")

    # Pre-generate iteration-specific seeds for cross-platform reproducibility
    # Each iteration gets its own deterministic seed regardless of which worker runs it
    iteration_seeds <- NULL
    if (is.null(set_seed)) {
        cli::cli_warn("No seed was set. Results may not be reproducible.")
    } else {
        # Save current RNG state to restore later
        old_rng_kind <- RNGkind()
        on.exit(RNGkind(old_rng_kind[1], old_rng_kind[2], old_rng_kind[3]), add = TRUE)
        
        # Set explicit RNG algorithm for cross-platform reproducibility
        # Mersenne-Twister is the most widely used and well-tested algorithm
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        set.seed(set_seed)
        
        # Pre-generate seeds for each iteration
        # This ensures each iteration uses the same seed regardless of worker assignment
        iteration_seeds <- sample.int(.Machine$integer.max, num_iter)
    }

    # Check object class
    if (!inherits(physeq, "phyloseq")) {
        stop("Input must be a phyloseq object, not a data.frame")
    }

    dataframe <- as.data.frame(
        as.matrix(t(otu_table(physeq, taxa_are_rows = TRUE)))
    )

    # Parallel setup
    threads <- min(threads, availableCores())
    cl <- makeCluster(threads)
    on.exit(stopCluster(cl), add = TRUE)

    # Export needed objects/packages to workers
    # Note: iteration_seeds may be NULL when set_seed is NULL, which is handled in worker logic
    clusterExport(
        cl,
        varlist = c("dataframe", "depth_level", "iteration_seeds"),
        envir = environment()
    )

    # Run rarefactions in parallel with per-iteration seeds
    com_iter <- parLapply(cl, 1:num_iter, function(i) {
        # Set iteration-specific seed for cross-platform reproducibility
        # Each iteration uses its pre-assigned seed regardless of which worker runs it
        if (!is.null(iteration_seeds)) {
            RNGkind("Mersenne-Twister", "Inversion", "Rejection")
            set.seed(iteration_seeds[i])
        }
        df <- rrarefy(dataframe, sample = depth_level)
        df <- as.data.frame(df)
        rownames_to_column(df, "sample_id")
    })

    # Aggregate results
    mean_data <- bind_rows(com_iter) |>
        group_by(sample_id) |>
        summarise(across(everything(), mean), .groups = "drop") |>
        filter(rowSums(across(where(is.numeric))) >= depth_level) |>
        column_to_rownames("sample_id")

    # Remove ASVs/OTUs with zero total abundance using base R
    mean_data <- mean_data[, colSums(mean_data) > 0]

    return(mean_data)
}
