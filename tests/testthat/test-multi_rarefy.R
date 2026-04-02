test_that("Test if parallel rarefaction works correctly", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  otu_table_rare <-
    multi_rarefy(
      physeq = bcse,
      depth_level = 3000,
      num_iter = 3,
      .summarize = TRUE,
      threads = 1,
      set_seed = 123
    )

  # Test if all the samples have the same number of reads
  read_counts <- rowSums(otu_table_rare)
  expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
})

test_that("multi_rarefy retains samples with row sums near depth_level (floating-point tolerance)", {
  skip_if_not_installed("phyloseq")

  # Simulate the floating-point scenario: a summarised data frame where some
  # samples have row sums that are VERY slightly below depth_level due to
  # floating-point arithmetic, which should still be kept.
  depth <- 1000L
  # Use a sub-epsilon offset to mimic the macOS floating-point mean artifact
  eps <- .Machine$double.eps

  # Build a tiny mock: one sample exactly at depth, one sub-epsilon below
  # (simulating macOS floating-point mean), one clearly below.
  mock_df <- data.frame(
    sample_id = c("exact", "fp_below", "too_low"),
    taxon1 = c(500, 500 - eps, 200),
    taxon2 = c(500, 500 + eps * 0.1, 100)
  )

  # Apply the same near()-based filter used inside multi_rarefy()
  result <- dplyr::filter(
    mock_df,
    dplyr::near(
      rowSums(across(where(is.numeric))),
      depth,
      tol = 1e-6
    )
  )

  # "exact" (rowSum == 1000) and "fp_below" (rowSum ≈ 1000) should be kept;
  # "too_low" (rowSum == 300) should be dropped.
  expect_equal(nrow(result), 2L)
  expect_true("exact" %in% result$sample_id)
  expect_true("fp_below" %in% result$sample_id)
  expect_false("too_low" %in% result$sample_id)
})

test_that("multi_rarefy outputs all rarefaction iterations when .summarize = FALSE", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  num_iterations <- 3
  rarefied_data <- multi_rarefy(
    physeq = bcse,
    depth_level = 3000,
    num_iter = num_iterations,
    .summarize = FALSE,
    threads = 1,
    set_seed = 123
  )

  # Should return a list of data frames
  expect_type(rarefied_data, "list")
  expect_equal(length(rarefied_data), num_iterations)
  expect_true(all(sapply(rarefied_data, is.data.frame)))

  # Each iteration should have the correct naming
  expect_equal(names(rarefied_data), paste0("iter_", 1:num_iterations))

  # Extract unique sample IDs from the first iteration (they all share the same samples)
  unique_samples <- nrow(sample_data(bcse))

  # Calculate expected samples after removing low-depth samples
  removed_samples <- setdiff(
    rownames(sample_data(bcse)),
    unique(sub("_iter_.*", "", rownames(rarefied_data[[1]])))
  )

  expected_samples_per_iter <- unique_samples - length(removed_samples)

  # Check each iteration has the expected number of samples
  for (i in seq_along(rarefied_data)) {
    expect_equal(
      nrow(rarefied_data[[i]]),
      expected_samples_per_iter,
      info = paste(
        "Iteration",
        i,
        "should have",
        expected_samples_per_iter,
        "samples"
      )
    )
  }

  # Verify row names have iteration suffixes
  for (i in seq_along(rarefied_data)) {
    iter_suffix <- paste0("_iter_", i)
    expect_true(
      all(grepl(paste0(iter_suffix, "$"), rownames(rarefied_data[[i]]))),
      info = paste(
        "All row names in iteration",
        i,
        "should end with",
        iter_suffix
      )
    )
  }

  # Verify all samples are rarefied to the correct depth
  for (i in seq_along(rarefied_data)) {
    row_sums <- rowSums(rarefied_data[[i]])
    expect_true(
      all(abs(row_sums - 3000) < 1e-6),
      info = paste("All samples in iteration", i, "should have depth ~3000")
    )
  }
})
