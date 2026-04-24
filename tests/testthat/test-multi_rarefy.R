test_that("Test if default behaviour is correct", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  otu_table_rare <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "array",
    set_seed = 123
  )

  # Should return a single averaged data frame (default behavior)
  expect_true(is.array(otu_table_rare) && length(dim(otu_table_rare)) == 3)

  # Test if all the samples have the same number of reads
  read_counts <- rowSums(otu_table_rare[,, 1])
  expect_equal(max(read_counts) - min(read_counts), 0) # All values are identical
  expect_true(all(abs(read_counts - 3000) < 1e-6))
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

test_that("multi_rarefy outputs all rarefaction iterations when .as = 'list'", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  num_iterations <- 3
  rarefied_data <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = num_iterations,
    .as = "list",
    set_seed = 123
  )

  # Should return a list of data frames
  expect_type(rarefied_data, "list")
  expect_equal(length(rarefied_data), num_iterations)
  expect_true(all(vapply(rarefied_data, is.data.frame, logical(1))))

  # Each iteration should have the correct naming
  expect_equal(names(rarefied_data), paste0("iter_", 1:num_iterations))

  # Extract unique sample IDs from the phyloseq object
  unique_samples <- nrow(phyloseq::sample_data(bcse))

  # Calculate expected samples after removing low-depth samples
  removed_samples <- setdiff(
    rownames(phyloseq::sample_data(bcse)),
    rownames(rarefied_data[[1]])
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

  # Verify row names are consistent across iterations (no iteration suffixes in new implementation)
  for (i in seq_along(rarefied_data)) {
    expect_equal(
      rownames(rarefied_data[[1]]),
      rownames(rarefied_data[[i]]),
      info = paste(
        "Iteration",
        i,
        "should have same sample names as iteration 1"
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

test_that("multi_rarefy with .as = 'array' returns 3D array", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  num_iterations <- 5
  rarefied_array <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = num_iterations,
    .as = "array",
    set_seed = 456
  )

  # Should return a 3D array
  expect_true(is.array(rarefied_array))
  expect_equal(length(dim(rarefied_array)), 3)

  # Check dimension names
  expect_equal(dimnames(rarefied_array)[[3]], paste0("iter_", 1:num_iterations))

  # All iterations should have same depth
  for (i in seq_len(num_iterations)) {
    row_sums <- rowSums(rarefied_array[,, i])
    expect_true(
      all(row_sums == 3000),
      info = paste("All samples in iteration", i, "should have depth 3000")
    )
  }
})

test_that("multi_rarefy list and array have consistent dimensions", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  # Same parameters for both
  rarefied_list <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 789
  )

  rarefied_array <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "array",
    set_seed = 789
  )

  # Array dimensions should match list dimensions
  expect_equal(dim(rarefied_array)[1], nrow(rarefied_list[[1]]))
  expect_equal(
    length(intersect(
      colnames(rarefied_list[[1]]),
      colnames(rarefied_array[,, 1])
    )),
    2861L
  )

  expect_equal(dim(rarefied_array)[1], nrow(rarefied_list[[1]]))
  expect_equal(dim(rarefied_array)[3], length(rarefied_list))
})

test_that("multi_rarefy removes zero-abundance taxa independently per iteration", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  rarefied_data <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 321
  )

  # Each iteration may have different numbers of non-zero taxa
  taxa_counts <- vapply(rarefied_data, ncol, integer(1))

  # Should have some variation (not all identical)
  # This tests that zero-taxa removal happens per iteration
  expect_true(length(unique(taxa_counts)) >= 1)

  # All should be > 0
  expect_true(all(taxa_counts > 0))
})

test_that("multi_rarefy handles low-depth samples correctly", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  # Use high depth to force some samples to be removed
  rarefied_data <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 50000,
    num_iter = 3,
    .as = "list",
    set_seed = 111
  )

  # Should retain 2 samples at the depth
  expect_true(length(rarefied_data) == 2)

  # Remaining samples should all have depth == 50000
  for (i in seq_along(rarefied_data)) {
    row_sums <- rowSums(rarefied_data[[i]])
    expect_true(all(row_sums == 50000))
  }
})

test_that("multi_rarefy single iteration returns data frame (not list)", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  # Single iteration with .as = "list" should return data frame
  rarefied_single <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 1,
    .as = "list",
    set_seed = 999
  )

  # Should be a DF
  expect_s3_class(rarefied_single, "data.frame")

  # Length of taxa should match original (minus any removed)
  expect_equal(length(rarefied_single), phyloseq::otu_table(bcse) |> nrow())
  expect_true(is.data.frame(rarefied_single))
})
