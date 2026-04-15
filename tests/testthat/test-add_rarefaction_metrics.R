# tests/testthat/test-add_rarefaction_metrics.R
test_that("add_rarefaction_metrics works with phyloseq object", {
  skip_if_not_installed("phyloseq")

  # Create simple phyloseq object with 10 samples
  otu_mat <- matrix(sample(1:100, 50, replace = TRUE), nrow = 5, ncol = 10)
  rownames(otu_mat) <- paste0("OTU", 1:5)
  colnames(otu_mat) <- paste0("Sample", 1:10)

  sample_df <- data.frame(
    treatment = rep(c("A", "B"), 5),
    row.names = paste0("Sample", 1:10)
  )

  physeq_obj <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_df)
  )

  # Run function
  result <- add_rarefaction_metrics(physeq_obj)

  # Check it's still phyloseq
  expect_s4_class(result, "phyloseq")

  # Check new columns added
  sample_data <- as.data.frame(phyloseq::sample_data(result))
  expect_true("read_num" %in% colnames(sample_data))
  expect_true("singleton_num" %in% colnames(sample_data))
  expect_true("goods_cov" %in% colnames(sample_data))
  expect_true("outlier" %in% colnames(sample_data))

  # Check we still have 10 samples
  expect_equal(nrow(sample_data), 10)
})

test_that("add_rarefaction_metrics works with a data.frame input", {
  otu_df <- as.data.frame(
    as(phyloseq::otu_table(bcse), "matrix")
  )

  result <- add_rarefaction_metrics(data = otu_df)

  expect_s3_class(result, "data.frame")
  expect_true(all(
    c("read_num", "singleton_num", "goods_cov", "outlier") %in% names(result)
  ))
  expect_equal(nrow(result), nrow(otu_df))
})

test_that("add_rarefaction_metrics data.frame output has correct metric values", {
  otu_df <- as.data.frame(
    as(phyloseq::otu_table(bcse), "matrix")
  )

  result <- add_rarefaction_metrics(data = otu_df)

  expect_true(all(result$read_num > 0))
  expect_true(all(result$singleton_num >= 0))
  expect_true(all(result$goods_cov >= 0 & result$goods_cov <= 100))
})

test_that("add_rarefaction_metrics rejects invalid input", {
  expect_error(
    add_rarefaction_metrics(data = matrix(1:10, nrow = 2)),
    "phyloseq object or a data.frame"
  )
})
