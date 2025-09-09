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
#' @param set_seed An optional integer to set the random seed for reproducibility (default = 123).
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
#' data(GlobalPatterns, package = "phyloseq")
#' test_otutable_rare <- parallel_rarefy(
#'   physeq = GlobalPatterns,
#'   depth_level = 500,
#'   num_iter = 10,
#'   set_seed = 123
#' )
#' 
#' rowSums(test_otutable_rare)
#'
#' @export
parallel_rarefy <- function(physeq,
                         depth_level,
                         num_iter = 100,
                         threads = get_available_cores(),
                         set_seed = NULL) {
    
    cat("\nSeed used:", set_seed, "\n")
    if (!is.null(set_seed)) set.seed(set_seed)
    
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
    
    return(mean_data)
}