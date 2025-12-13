# tests/testthat/test-plot_core_distribution.R
test_that("plot_core_distribution returns a ggplot and builds", {
    testthat::skip_if_not_installed("ggplot2")

    # Create a matrix of random counts
    otu_table <- base::as.data.frame(
        matrix(
            rpois(6 * 15, lambda = 20) * (runif(6 * 15) > 0.4),
            nrow = 15,
            ncol = 6
        )
    )

    rownames(otu_table) <- paste0("OTU", 1:15)
    colnames(otu_table) <- paste0("Sample", 1:6)

    metadata <- base::data.frame(
        sample_id = colnames(otu_table),
        time = factor(c(rep(2016, 3), rep(2017, 3)))
    )

    increase_core <- c("OTU1", "OTU2", "OTU3", "OTU5")

    # Create the dataframe
    otu_ranked <- data.frame(
        otu = rownames(otu_table),
        abundance = rowSums(otu_table)
    ) %>%
        dplyr::arrange(desc(abundance)) %>%
        mutate(rank = row_number())

    core_result_example <- list(
        otu_table = otu_table,
        metadata = metadata,
        increase_core = increase_core,
        otu_ranked = otu_ranked
    )

    p <- plot_core_distribution(
        core_result = core_result_example,
        core_set = "increase",
        group_var = "time"
    )

    testthat::expect_s3_class(p, "ggplot")
    testthat::expect_silent(ggplot2::ggplot_build(p))
})
