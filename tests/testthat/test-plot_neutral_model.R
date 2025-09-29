# tests/testthat/test-plot_neutral_model.R
test_that("plot_neutral_model returns a ggplot and builds", {
    testthat::skip_if_not_installed("ggplot2")
    
    fit_result <- list(
        model_prediction = data.frame(
            otu_rel    = c(1e-4, 5e-4, 1e-3, 2e-3, 1e-2),
            otu_occ    = c(0.10, 0.20, 0.40, 0.60, 0.80),
            membership = c(NA, "core", "core", NA, "core"),
            fit_class  = c("As predicted","Above prediction",
                           "Below prediction","As predicted","As predicted"),
            p          = c(1e-4, 5e-4, 1e-3, 2e-3, 1e-2),
            freq.pred  = c(0.12, 0.25, 0.38, 0.62, 0.79),
            pred.lwr   = c(0.05, 0.15, 0.30, 0.50, 0.70),
            pred.upr   = c(0.20, 0.35, 0.50, 0.75, 0.90)
        ),
        goodness_of_fit = data.frame(Rsqr = 0.83, m = 0.159)
    )
    
    p <- plot_neutral_model(fit_result)
    testthat::expect_s3_class(p, "ggplot")
    testthat::expect_silent(ggplot2::ggplot_build(p))
})