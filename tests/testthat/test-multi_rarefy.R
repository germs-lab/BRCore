test_that("rarefaction worked correctly", {
    data(esophagus, package = "phyloseq")
    
    otu_table_rare <<-
        BRCore::multi_rarefy(physeq = test_phyloseq2,
                             depth_level = 500,
                             num_iter = 10)
    
    # Test if all the samples have the same number of reads
    expect_true(all(rowSums(otu_table_rare) == rowSums(otu_table_rare)[1]))
})