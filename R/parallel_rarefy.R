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
#' @import phyloseq
#' @import vegan
#' @import dplyr
#' @import tibble
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom dplyr %>% 
#'
#' @examples
#' \donttest{
#' library(phyloseq)
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Example rarefaction (single iteration, single core to keep examples fast)
#' otu_table_rare <-
#'     parallel_rarefy(
#'         physeq = GlobalPatterns,
#'         depth_level = 200,
#'         num_iter = 3,
#'         threads = 1,
#'         set_seed = 123
#'     )
#'
#' rowSums(otu_table_rare)
#' }
#'
#' @export
parallel_rarefy <- function(physeq,
                         depth_level,
                         num_iter = 100,
                         threads = get_available_cores(),
                         set_seed = NULL) {
    
    cli::cli_text("\nSeed used: {set_seed}\n")
    
    if (is.null(set_seed)) {
        cli::cli_warn("No seed was set. Results may not be reproducible.")
    } else {
        set.seed(set_seed)
    }
  
    
    # Check object class
    if (!inherits(physeq, "phyloseq")) {
        stop("Input must be a phyloseq object, not a data.frame")
    }
    
    dataframe <- as.data.frame(
        as.matrix(t(phyloseq::otu_table(physeq, taxa_are_rows = TRUE)))
    )
    
    # Parallel setup
    threads <- min(threads, parallelly::availableCores())
    cl <- parallel::makeCluster(threads)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Export needed objects/packages to workers
    parallel::clusterExport(
        cl,
        varlist = c("dataframe", "depth_level"),
        envir = environment()
    )
    
    # Run rarefactions in parallel
    com_iter <- parallel::parLapply(cl, 1:num_iter, function(i) {
        set.seed(set_seed + i) # each worker gets a different but reproducible seed
        df <- vegan::rrarefy(dataframe, sample = depth_level)
        df <- as.data.frame(df)
        tibble::rownames_to_column(df, "sample_id")
    })
    
    # Aggregate results
    mean_data <- dplyr::bind_rows(com_iter) %>%
        dplyr::group_by(sample_id) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
        dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) >= depth_level) %>%
        tibble::column_to_rownames("sample_id")
    
    # Remove ASVs/OTUs with zero total abundance using base R
    mean_data <- mean_data[, colSums(mean_data) > 0]
    
    return(mean_data)
}