test_that("phyloseq object correctly created", {
    # Load expected
    #load(here::here("tests/testthat/testdata/test_phyloseq.rda"))
    load(testthat::test_path("testdata", "test_phyloseq.rda"))
    
    otu_table_rare <-
        multi_rarefy(
            physeq = test_phyloseq,
            depth_level = 500,
            num_iter = 10,
            set_seed = 1234
        )
    
    test_phyloseq_rare <- do_phyloseq(test_phyloseq, otu_table_rare)
    read_counts <- sample_sums(test_phyloseq_rare)
    
    
    expect_true(class(test_phyloseq_rare) == "phyloseq")
    
    expect_true(is.numeric(read_counts))
    expect_true(length(read_counts) > 0)
    expect_true(all(read_counts > 0))
    
    # Test if all the samples have the same number of reads
    expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
    expect_equal(median(read_counts), mean(read_counts))
    
    
})
