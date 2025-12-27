# tests/testthat/test-add_rarefaction_metrics.R
test_that("add_rarefaction_metrics works with phyloseq object", {
    skip_if_not_installed("phyloseq")
    
    # Create simple phyloseq object with 10 samples
    otu_mat <- matrix(sample(1:100, 50, replace = TRUE), nrow = 5, ncol = 10)
    rownames(otu_mat) <- paste0("OTU", 1:5)
    colnames(otu_mat) <- paste0("Sample", 1:10)
    
    sample_df <- data.frame(
        treatment = rep(c("A", "B"), 5),
        row.names = paste0("Sample", 1:10)
    )
    
    physeq <- phyloseq::phyloseq(
        phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
        phyloseq::sample_data(sample_df)
    )
    
    # Run function
    result <- add_rarefaction_metrics(physeq)
    
    # Check it's still phyloseq
    expect_s4_class(result, "phyloseq")
    
    # Check new columns added
    sample_data <- as.data.frame(phyloseq::sample_data(result))
    expect_true("read_num" %in% colnames(sample_data))
    expect_true("singlton_num" %in% colnames(sample_data))
    expect_true("goods_cov" %in% colnames(sample_data))
    expect_true("outlier" %in% colnames(sample_data))
    
    # Check we still have 10 samples
    expect_equal(nrow(sample_data), 10)
})