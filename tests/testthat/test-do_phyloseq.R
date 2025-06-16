test_that("phyloseq object correctly created", {
    # Load expected
    load(here::here("tests/data/test_phyloseq.rda"))
    
    
    otu_table_rare <-
        multi_rarefy(
            physeq = test_phyloseq,
            depth_level = 500,
            num_iter = 10,
            set_seed = 1234
        )
    
    test_phyloseq_rare <- do_phyloseq(test_phyloseq, otu_table_rare)
    read_counts <- sample_sums(test_phyloseq_rare)
    
    # Test if all the samples have the same number of reads
    expect_true(all(read_counts == read_counts[1]))
    # Test if the created phyloseq object is class phyloseq
    expect_true(class(test_phyloseq_rare) == "phyloseq")
})