test_that("rarefaction worked correctly", {
    # Load expected
    #load(here::here("tests/testthat/testdata/test_phyloseq.rda"))
    load(testthat::test_path("testdata", "test_phyloseq.rda"))
    
    otu_table_rare <-
        parallel_rarefy(
            physeq = test_phyloseq,
            depth_level = 200,
            num_iter = 3,
            threads = 1,
            set_seed = 123
        )
    
    # Test if all the samples have the same number of reads
    read_counts <- rowSums(otu_table_rare)
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
    
})
