#' Perform Multiple Rarefaction on a Phyloseq Object
#'
#' This function performs multiple rarefaction on a `phyloseq` object by randomly sub-sampling OTUs/ASVs within samples without replacement. The process is repeated for a specified number of iterations, and the results are averaged. Samples with fewer OTUs/ASVs than the specified `depth_level` are discarded.
#'
#' @param physeq A `phyloseq` object containing an OTU/ASV table.
#' @param depth_level An integer specifying the sequencing depth (number of OTUs/ASVs) to which samples should be rarefied.
#' @param num_iter An integer specifying the number of iterations to perform for rarefaction.
#'
#' @return A data frame with samples as rows and taxa as columns. The values represent the average abundance across all iterations.
#'
#' @importFrom phyloseq otu_table
#' @importFrom vegan rrarefy
#' @import dplyr
#' @import tibble
#' @export
#'
#' @examples
#' # Load example data
#' data(GlobalPatterns, package = "phyloseq")
#'
#' # Perform multiple rarefaction
#' rarefied_data <- multi_rarefy(GlobalPatterns, depth_level = 1000, num_iter = 10)
#' head(rarefied_data)
multi_rarefy <- function(physeq, depth_level, num_iter) {

  dataframe <-
    as.data.frame(as.matrix(t(otu_table(physeq, taxa_are_rows = TRUE))))
  
  com_iter <- vector(mode = "list", length =  num_iter)
  
  for (i in seq_along(com_iter)) {
    com_iter[[i]] <- as.data.frame(vegan::rrarefy(dataframe, sample = depth_level)) %>% rownames_to_column("SampleID")
  }
  
  mean_data <- do.call(rbind, com_iter)
  mean_data <- mean_data %>%
    group_by(SampleID) %>%
    summarise(across(everything(), mean)) %>%
    filter(rowSums(across(where(is.numeric))) >= depth_level)
  
  print(mean_data %>% as_tibble())
  
  return(mean_data)
}






