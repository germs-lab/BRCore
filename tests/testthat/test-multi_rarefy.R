test_that("rarefaction worked correctly", {
    # Load expected
    load(here::here("tests/data/test_phyloseq.rda"))
    
    
    otu_table_rare <-
        multi_rarefy(
            physeq = test_phyloseq,
            depth_level = 500,
            num_iter = 10,
            set_seed = 1234
        )
    
    # Test if all the samples have the same number of reads
    expect_true(all(rowSums(otu_table_rare) == 500))
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
    
})
