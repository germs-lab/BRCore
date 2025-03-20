# Test suite
test_that("calculate_nmds() works with esophagus dataset", {
  # Load the esophagus dataset
  data(esophagus, package = "phyloseq")
  
  # Add sample metadata to the esophagus dataset (required for calculate_nmds)
  sample_data <- data.frame(
    Sample = phyloseq::sample_names(esophagus),
    Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE), # Random groups
    row.names = phyloseq::sample_names(esophagus)
  )
  esophagus_with_metadata <- phyloseq::merge_phyloseq(esophagus, sample_data(sample_data))
  
  # Load expected results
  load(here::here("tests/testthat/expected_nmds.rda"))

  # Extract the OTU table from the esophagus dataset
  esophagus_asv_matrix <- as.matrix(otu_table(esophagus_with_metadata))

  # Run the function
  results <- calculate_nmds(esophagus_asv_matrix, esophagus_with_metadata,
    ncores = 1,
    k = 2,
    trymax = 100
  )

  # Test 1: Check output structure
  expect_named(results, c("nmds_scores", "nmds_df", "ordi"))

  # Test 2: Check nmds_scores
  expect_s3_class(results$nmds_scores, "data.frame")
  expect_named(results$nmds_scores, c("unique_id", "NMDS1", "NMDS2"))
  expect_equal(nrow(results$nmds_scores), nrow(esophagus_with_metadata@otu_table))

  # Test 3: Check nmds_df
  expect_s3_class(results$nmds_df, "data.frame")
  expect_true(all(c("unique_id", "NMDS1", "NMDS2", "Sample", "Group") %in% names(results$nmds_df)))

  # Test 4: Check ordi
  expect_s3_class(results$ordi, "metaMDS")
  expect_equal(results$ordi$ndim, 2) # Check number of dimensions

  # Test 5: Compare NMDS scores with expected values (using waldo)
  expected_nmds_scores <- expected_nmds$nmds_scores
  compare(results$nmds_scores, expected_nmds_scores)

  # Test 6: Compare merged metadata with expected values (using waldo)
  expected_nmds_df <- expected_nmds$nmds_df
  compare(results$nmds_df, expected_nmds_df)

  # Test 7: Check for NA/NaN in NMDS scores
  expect_false(any(is.na(results$nmds_scores)))

  # Test 8: Check for convergence
  expect_true(ifelse(results$ordi$converged > 0, TRUE, FALSE))
})
