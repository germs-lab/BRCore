test_that("GlobalPatterns pipeline -> fit_neutral_model basic structure", {
  skip_on_cran()
  skip_if_not(exists("identify_core"), "identify_core() not found")
  skip_if_not(exists("fit_neutral_model"), "fit_neutral_model() not found")

  test_that("switchgrass dataset loads", {
    utils::data("switchgrass", package = "BRCore", envir = environment())
    expect_true(exists("switchgrass"))
    # e.g., check structure
    expect_s4_class(switchgrass, "phyloseq")
  })

  switchgrass_core <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 123
  )

  switchgrass_core_fit <- fit_neutral_model(
    otu_table = switchgrass_core$otu_table,
    core_set = switchgrass_core$increase_core,
    abundance_occupancy = switchgrass_core$abundance_occupancy
  )

  expect_type(switchgrass_core_fit, "list")
  expect_true(all(
    c("model_prediction", "goodness_of_fit") %in%
      names(switchgrass_core_fit)
  ))

  if (exists("plot_neutral_model")) {
    p <- plot_neutral_model(switchgrass_core_fit)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("fit_neutral_model errors on invalid otu_table type", {
  expect_error(
    fit_neutral_model(
      otu_table = list(a = 1, b = 2),
      core_set = "ASV1",
      abundance_occupancy = data.frame(otu = "ASV1")
    ),
    "matrix or data.frame"
  )
})

test_that("fit_neutral_model errors on NULL core_set", {
  otu <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(
      c("S1", "S2", "S3"),
      c("ASV1", "ASV2", "ASV3", "ASV4")
    )
  )
  expect_error(
    fit_neutral_model(
      otu_table = otu,
      core_set = NULL,
      abundance_occupancy = data.frame(otu = c("ASV1", "ASV2"))
    ),
    "non-empty character vector"
  )
})

test_that("fit_neutral_model errors on empty core_set", {
  otu <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(
      c("S1", "S2", "S3"),
      c("ASV1", "ASV2", "ASV3", "ASV4")
    )
  )
  expect_error(
    fit_neutral_model(
      otu_table = otu,
      core_set = character(0),
      abundance_occupancy = data.frame(otu = c("ASV1", "ASV2"))
    ),
    "non-empty character vector"
  )
})

test_that("fit_neutral_model errors when abundance_occupancy lacks 'otu' column", {
  otu <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(
      c("S1", "S2", "S3"),
      c("ASV1", "ASV2", "ASV3", "ASV4")
    )
  )
  expect_error(
    fit_neutral_model(
      otu_table = otu,
      core_set = c("ASV1"),
      abundance_occupancy = data.frame(taxon = c("ASV1", "ASV2"))
    ),
    "otu"
  )
})
