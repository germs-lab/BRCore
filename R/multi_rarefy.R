# The functions below runs a multiple rarefaction. 
# Multiple rarefaction performs random sub-sampling without replacement of OTU/ASVs 
# within samples, repeated for a specified number of iterations, and then 
# averages the results. Samples with less OTU/ASVs than the specified depth_level
# are not rarefied and are discarded from the final object.

# The function output a data frame with samples as rows and taxa as columns.

multi_rarefy <- function(physeq, depth_level, num_iter){
  require(tidyverse)
  require(vegan)
  
  dataframe <- 
    as.data.frame(as.matrix(t(physeq@otu_table)))
  
  com_iter <- vector(mode = "list", length =  num_iter)
  
  for (i in seq_along(com_iter)) {
    com_iter[[i]] <- as.data.frame(
      vegan::rrarefy(dataframe, sample = depth_level)
    ) %>% rownames_to_column("SampleID")
  }
  
  mean_data <- do.call(rbind, com_iter)
  mean_data <- mean_data %>% 
    group_by(SampleID) %>% 
    summarise(across(everything(), mean)) %>% 
    filter(rowSums(across(where(is.numeric))) >= depth_level)
  
  print(mean_data %>% as_tibble())
  
  return(mean_data)
}






