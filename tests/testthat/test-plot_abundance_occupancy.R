# tests/testthat/test-plot_abundance_occupancy.R
test_that("plot_abundance_occupancy returns a ggplot and builds", {
    testthat::skip_if_not_installed("ggplot2")
    
    abundance_occupancy <- data.frame(
        otu = paste0("OTU", 1:10),
        otu_occ = round(runif(10, min = 0.05, max = 1), 8),
        otu_rel = round(runif(10, min = 0.0002, max = 0.05), 8),
        membership = sample(c("Core", "Not core"), size = 10, replace = TRUE, prob = c(0.4, 0.6))
    )
    
    increase_core <- c("OTU1", "OTU4", "OTU7")
    
    core_result_example <- list(
        abundance_occupancy = abundance_occupancy,
        increase_core = "incrase"
    )
    
    p <- plot_abundance_occupancy(core_result = core_result_example,
                                  core_set = "increase")
    
    testthat::expect_s3_class(p, "ggplot")
    testthat::expect_silent(ggplot2::ggplot_build(p))
})