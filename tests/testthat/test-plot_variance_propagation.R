# tests/testthat/test-plot_variance_propagation.R

test_that("plot_variance_propagation rejects non-phyloseq input", {
  skip_if_not_installed("phyloseq")

  expect_error(
    plot_variance_propagation(
      physeq_obj = data.frame(a = 1),
      rarefied = list(),
      q = 0,
      group_var = "x",
      group_color = "y"
    ),
    "phyloseq"
  )
})

# List input ----

test_that("plot_variance_propagation returns ggplot with list input", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 7643
  )

  p <- plot_variance_propagation(
    physeq_obj = bcse,
    rarefied = rarefied,
    q = 0,
    group_var = "Crop",
    group_color = "Plot"
  )

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plot works with q = 1 (Shannon)", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  p <- plot_variance_propagation(
    physeq_obj = bcse,
    rarefied = rarefied,
    q = 1,
    group_var = "Crop",
    group_color = "Plot"
  )

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plot works with q = 2 (Simpson)", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  p <- plot_variance_propagation(
    physeq_obj = bcse,
    rarefied = rarefied,
    q = 2,
    group_var = "Crop",
    group_color = "Plot"
  )

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

# Array input ----

test_that("plot_variance_propagation returns ggplot with array input", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "array",
    set_seed = 42
  )

  expect_true(is.array(rarefied))
  expect_equal(length(dim(rarefied)), 3)

  p <- plot_variance_propagation(
    physeq_obj = bcse,
    rarefied = rarefied,
    q = 1,
    group_var = "Crop",
    group_color = "Plot"
  )

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

# Plot structure ----

test_that("plot has correct facets (Raw vs Rarefied)", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  p <- plot_variance_propagation(
    physeq_obj = bcse,
    rarefied = rarefied,
    q = 0,
    group_var = "Crop",
    group_color = "Plot"
  )

  build <- ggplot2::ggplot_build(p)
  plot_data <- build$data[[1]]

  expect_true("PANEL" %in% names(plot_data))
  expect_equal(length(unique(plot_data$PANEL)), 2)
})

test_that("plot y-axis label reflects q value", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  for (q_val in c(0, 1, 2)) {
    p <- plot_variance_propagation(
      physeq_obj = bcse,
      rarefied = rarefied,
      q = q_val,
      group_var = "Crop",
      group_color = "Plot"
    )

    expect_equal(p$labels$y, paste0("Hill number (q = ", q_val, ")"))
  }
})

# Sample matching ----

test_that("plot handles rarefied subset of samples", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("ggplot2")

  data("bcse", package = "BRCore")

  bcse_sub <- phyloseq::subset_samples(bcse, Crop == "Switchgrass")

  rarefied <- multi_rarefy(
    physeq_obj = bcse_sub,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  p <- plot_variance_propagation(
    physeq_obj = bcse_sub,
    rarefied = rarefied,
    q = 1,
    group_var = "Crop",
    group_color = "Plot"
  )

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

# Vector q (currently errors) ----

test_that("plot_variance_propagation errors with vector q", {
  skip_if_not_installed("phyloseq")

  data("bcse", package = "BRCore")

  rarefied <- multi_rarefy(
    physeq_obj = bcse,
    depth_level = 3000,
    num_iter = 3,
    .as = "list",
    set_seed = 42
  )

  # .hill_number does not support length(q) > 1
  expect_error(
    plot_variance_propagation(
      physeq_obj = bcse,
      rarefied = rarefied,
      q = c(0, 1, 2),
      group_var = "Crop",
      group_color = "Plot"
    )
  )
})
