# tests/testthat/test-plot_rarefaction_metrics.R
test_that("plot_rarefaction_metrics returns a ggplot and builds", {
    testthat::skip_if_not_installed("ggplot2")
    
    rare_metrics <- data.frame(
        otu = paste0("sample", 1:10),
        read_num = c(25, 17, 6002, 6428, 8854, 8281, 6707, 7133, 7559, 8985),
        singlton_num = c(117, 137, 26, 37, 38, 40, 23, 32, 34, 55),
        goods_cov = c(0.719, 0.762, 0.988, 0.987, 0.826, 0.972, 0.909, 0.875, 0.885, 0.934),
        outlier = c(25, 17, NA, NA, NA, NA, NA, NA, NA, NA)
    )
    
    p <- plot_rarefaction_metrics(data = rare_metrics)
    
    testthat::expect_s3_class(p, "ggplot")
    testthat::expect_silent(ggplot2::ggplot_build(p))
})