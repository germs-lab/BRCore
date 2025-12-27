test_that("phyloseq object correctly created", {
    skip_if_not_installed("phyloseq")

    data("bcse", package = "BRCore")

    otu_table_rare <-
        multi_rarefy(
            physeq = bcse,
            depth_level = 200,
            num_iter = 3,
            threads = 1,
            set_seed = 123
        )

    test_bcse_rare <- update_otu_table(bcse, otu_table_rare)
    read_counts <- sample_sums(test_bcse_rare)

    expect_true(class(test_bcse_rare) == "phyloseq")

    expect_true(is.numeric(read_counts))
    expect_true(length(read_counts) > 0)
    expect_true(all(read_counts > 0))

    # Test if all the samples have the same number of reads
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
    expect_equal(median(read_counts), mean(read_counts))
})
