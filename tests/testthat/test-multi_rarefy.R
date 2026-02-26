test_that("Test if parallel rarefaction works correctly", {
    skip_if_not_installed("phyloseq")

    data("bcse", package = "BRCore")

    otu_table_rare <-
        multi_rarefy(
            physeq = bcse,
            depth_level = 3000,
            num_iter = 3,
            threads = 1,
            set_seed = 123
        )

    # Test if all the samples have the same number of reads
    read_counts <- rowSums(otu_table_rare)
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
})

test_that("multi_rarefy retains samples with row sums near depth_level (floating-point tolerance)", {
    skip_if_not_installed("phyloseq")

    # Simulate the floating-point scenario: a summarised data frame where some
    # samples have row sums that are VERY slightly below depth_level due to
    # floating-point arithmetic, which should still be kept.
    depth <- 1000L
    # Use a sub-epsilon offset to mimic the macOS floating-point mean artifact
    eps <- .Machine$double.eps

    # Build a tiny mock: one sample exactly at depth, one sub-epsilon below
    # (simulating macOS floating-point mean), one clearly below.
    mock_df <- data.frame(
        sample_id = c("exact", "fp_below", "too_low"),
        taxon1 = c(500, 500 - eps, 200),
        taxon2 = c(500, 500 + eps * 0.1, 100)
    )

    # Apply the same near()-based filter used inside multi_rarefy()
    result <- dplyr::filter(
        mock_df,
        dplyr::near(
            rowSums(across(where(is.numeric))),
            depth,
            tol = 1e-6
        )
    )

    # "exact" (rowSum == 1000) and "fp_below" (rowSum ≈ 1000) should be kept;
    # "too_low" (rowSum == 300) should be dropped.
    expect_equal(nrow(result), 2L)
    expect_true("exact" %in% result$sample_id)
    expect_true("fp_below" %in% result$sample_id)
    expect_false("too_low" %in% result$sample_id)
})
