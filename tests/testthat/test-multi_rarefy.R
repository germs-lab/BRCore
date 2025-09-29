test_that("Test if parallel rarefaction works correctly", {
    # Load data
    load(testthat::test_path("testdata", "test_phyloseq.rda"))
    
    otu_table_rare <-
        multi_rarefy(
            physeq = test_phyloseq,
            depth_level = 300,
            num_iter = 3,
            threads = 1,
            set_seed = 123
        )
    
    # Test if all the samples have the same number of reads
    read_counts <- rowSums(otu_table_rare)
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
    
})
