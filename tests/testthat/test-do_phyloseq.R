test_that("phyloseq object correctly created", {
    test_phyloseq3 <<- BRCore::do_phyloseq(test_phyloseq2, otu_table_rare)
    read_counts <- sample_sums(test_phyloseq3)
    # Test if all the samples have the same number of reads
    expect_true(all(read_counts == read_counts[1]))
    # Test if the created phyloseq object is class phyloseq
    expect_true(class(test_phyloseq3) == "phyloseq")
})