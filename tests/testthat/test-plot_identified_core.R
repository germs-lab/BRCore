# tests/testthat/test-plot_identified_core.R
test_that("plot_identified_core returns a ggplot and builds", {
    skip_if_not_installed("ggplot2")
    skip_if_not(
        exists("plot_identified_core"),
        "plot_identified_core() not found"
    )

    # minimal deterministic mock input
    n <- 80
    bray_curtis_ranked <- data.frame(
        rank = as.character(seq_len(n)), # function coerces this
        proportionBC = seq(0.01, 0.9, length.out = n)
    )
    elbow <- 10
    lastCall <- 30

    p <- plot_identified_core(
        bray_curtis_ranked,
        elbow = elbow,
        lastCall = lastCall
    )

    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p)) # renders without error
})
