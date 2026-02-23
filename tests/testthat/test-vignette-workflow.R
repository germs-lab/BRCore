test_that("1: Vignette test data structures are valid", {
  skip_on_cran()

  # Load reference data
  load(here("tests/test_sets/test_vignette_data.rda"))

  # Test: Rarefied OTU table has correct dimensions and sums
  expect_true(
    is.matrix(test_vignette_data$test_bcse_rarefied_otutable) ||
      is.data.frame(test_vignette_data$test_bcse_rarefied_otutable)
  )

  expect_equal(999.999999999996, 1000, tol = 1e-6) # Sanity check for floating-point precision

  row_sums <- rowSums(test_vignette_data$test_bcse_rarefied_otutable)
  expect_true(all(near(row_sums, 1000, tol = 1e-6)))

  # Test: Core result structure is complete
  expect_named(
    test_vignette_data$test_bcse_rare_core,
    c(
      "bray_curtis_ranked",
      "otu_ranked",
      "abundance_occupancy",
      "priority_var",
      "abundance_weight",
      "elbow",
      "bc_increase",
      "increase_value",
      "elbow_core",
      "increase_core",
      "otu_table",
      "metadata",
      "taxonomy"
    )
  )

  # Test: Core sets are non-empty
  expect_true(length(test_vignette_data$test_bcse_rare_core$elbow_core) > 0)
  expect_true(length(test_vignette_data$test_bcse_rare_core$increase_core) > 0)

  # Test: Increase core is larger than elbow core
  expect_true(
    length(test_vignette_data$test_bcse_rare_core$increase_core) >=
      length(test_vignette_data$test_bcse_rare_core$elbow_core)
  )

  # Test: All plots are ggplot objects
  plot_objects <- c(
    "test_rarefaction_plot",
    "test_bcse_identified_core",
    "test_plot_abund_occ_increase",
    "test_plot_abund_occ_elbow",
    "test_plot_core_dist_bar",
    "test_plot_core_dist_line",
    "test_plot_core_dist_heatmap",
    "test_plot_bcse_neutral_fit"
  )

  for (plot_name in plot_objects) {
    expect_s3_class(test_vignette_data[[plot_name]], "ggplot")
  }
})


test_that("2: Vignette workflow produces consistent results", {
  skip_on_cran()
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("tidyverse")
  skip_if_not_installed("viridis")

  # Load reference data from saved vignette run
  load(here("tests/test_sets/test_vignette_data.rda"))

  # Run the vignette workflow step by step
  library(phyloseq)
  library(tidyverse)
  library(viridis)

  # Step 1: Load bcse data
  data("bcse", package = "BRCore")

  # Step 2: Calculate rarefaction metrics
  bcse_metrics <- add_rarefaction_metrics(data = bcse)

  # Test: Rarefaction metrics match
  expect_s4_class(bcse_metrics, "phyloseq")
  expect_identical(
    sample_data(bcse_metrics),
    sample_data(test_vignette_data$test_bcse_metrics)
  )

  # Step 3: Plot rarefaction metrics
  rarefaction_plot <- plot_rarefaction_metrics(bcse_metrics)

  # Test: Plot structure matches
  expect_s3_class(rarefaction_plot, "ggplot")
  expect_identical(
    rarefaction_plot$labels,
    test_vignette_data$test_rarefaction_plot$labels
  )

  # Step 4: Multiple rarefaction
  # Use single thread on CI to ensure reproducibility
  n_threads <- if (nzchar(Sys.getenv("CI"))) {
    1
  } else if (.Platform$OS.type == "windows") {
    1
  } else {
    2
  }

  bcse_rarefied_otutable <- multi_rarefy(
    physeq = bcse,
    depth_level = 1000,
    num_iter = 100,
    threads = n_threads,
    set_seed = 7642
  )

  # Test: Rarefied OTU table properties match
  expect_equal(
    dim(bcse_rarefied_otutable),
    dim(test_vignette_data$test_bcse_rarefied_otutable)
  )
  # Test: All samples have row sums of 1000 (accounting for floating-point precision)
  expect_equal(
    rowSums(bcse_rarefied_otutable),
    setNames(
      rep(1000, nrow(bcse_rarefied_otutable)),
      rownames(bcse_rarefied_otutable)
    ),
    tolerance = 1e-6
  )
  expect_equal(
    bcse_rarefied_otutable,
    test_vignette_data$test_bcse_rarefied_otutable,
    tolerance = 1e-6
  )

  # Step 5: Update OTU table
  bcse_rare <- update_otu_table(
    physeq = bcse,
    otu_rare = bcse_rarefied_otutable
  )

  # Test: Updated phyloseq object matches
  expect_s4_class(bcse_rare, "phyloseq")
  expect_equal(ntaxa(bcse_rare), ntaxa(test_vignette_data$test_bcse_rare))
  expect_equal(nsamples(bcse_rare), nsamples(test_vignette_data$test_bcse_rare))
  expect_equal(
    sample_sums(bcse_rare),
    sample_sums(test_vignette_data$test_bcse_rare),
    tolerance = 1e-6
  )

  # Step 6: Identify core microbiome
  bcse_rare_core <- identify_core(
    physeq_obj = bcse_rare,
    priority_var = "Crop",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 2134
  )

  # Test: Core identification results match
  expect_s3_class(bcse_rare_core, "identify_core_result")
  expect_equal(
    bcse_rare_core$priority_var,
    test_vignette_data$test_bcse_rare_core$priority_var
  )
  expect_equal(
    bcse_rare_core$abundance_weight,
    test_vignette_data$test_bcse_rare_core$abundance_weight
  )
  expect_equal(
    bcse_rare_core$elbow,
    test_vignette_data$test_bcse_rare_core$elbow
  )
  expect_equal(
    bcse_rare_core$bc_increase,
    test_vignette_data$test_bcse_rare_core$bc_increase
  )
  expect_equal(
    bcse_rare_core$increase_value,
    test_vignette_data$test_bcse_rare_core$increase_value
  )

  # Test: Core sets match
  expect_equal(
    bcse_rare_core$elbow_core,
    test_vignette_data$test_bcse_rare_core$elbow_core
  )
  expect_equal(
    bcse_rare_core$increase_core,
    test_vignette_data$test_bcse_rare_core$increase_core
  )

  # Test: Bray-Curtis ranked data matches
  expect_equal(
    bcse_rare_core$bray_curtis_ranked,
    test_vignette_data$test_bcse_rare_core$bray_curtis_ranked,
    tolerance = 1e-6
  )

  # Step 7: Plot identified core
  bcse_identified_core <- plot_identified_core(
    bray_curtis_ranked = bcse_rare_core$bray_curtis_ranked,
    elbow = bcse_rare_core$elbow,
    lastCall = bcse_rare_core$bc_increase,
    increase_value = bcse_rare_core$increase_value
  )

  # Test: Core identification plot structure
  expect_s3_class(bcse_identified_core, "ggplot")
  expect_identical(
    bcse_identified_core$labels,
    test_vignette_data$test_bcse_identified_core$labels
  )

  # Step 8: Plot abundance-occupancy (increase)
  plot_abund_occ_increase <- plot_abundance_occupancy(
    core_result = bcse_rare_core,
    core_set = "increase"
  )

  # Test: Abundance-occupancy plot structure (increase)
  expect_s3_class(plot_abund_occ_increase, "ggplot")
  expect_identical(
    plot_abund_occ_increase$labels,
    test_vignette_data$test_plot_abund_occ_increase$labels
  )

  # Step 9: Plot abundance-occupancy (elbow)
  plot_abund_occ_elbow <- plot_abundance_occupancy(
    core_result = bcse_rare_core,
    core_set = "elbow"
  )

  # Test: Abundance-occupancy plot structure (elbow)
  expect_s3_class(plot_abund_occ_elbow, "ggplot")
  expect_identical(
    plot_abund_occ_elbow$labels,
    test_vignette_data$test_plot_abund_occ_elbow$labels
  )
})


test_that("3: Vignette core distribution plots are consistent", {
  skip_on_cran()
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("tidyverse")
  skip_if_not_installed("viridis")

  # Load reference data
  load(here("tests/test_sets/test_vignette_data.rda"))

  # Recreate necessary objects for plotting
  data("bcse", package = "BRCore")
  library(phyloseq)
  library(tidyverse)
  library(viridis)

  # Use single thread on CI to ensure reproducibility
  n_threads <- if (nzchar(Sys.getenv("CI"))) {
    1
  } else if (.Platform$OS.type == "windows") {
    1
  } else {
    2
  }

  bcse_rarefied_otutable <- multi_rarefy(
    physeq = bcse,
    depth_level = 1000,
    num_iter = 100,
    threads = n_threads,
    set_seed = 7642
  )

  bcse_rare <- update_otu_table(
    physeq = bcse,
    otu_rare = bcse_rarefied_otutable
  )

  bcse_rare_core <- identify_core(
    physeq_obj = bcse_rare,
    priority_var = "Crop",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 2134
  )

  # Test bar plot
  plot_core_dist_bar <- plot_core_distribution(
    core_result = bcse_rare_core,
    core_set = "increase",
    group_var = "Crop",
    plot_type = "bar"
  )

  expect_s3_class(plot_core_dist_bar, "ggplot")
  expect_identical(
    plot_core_dist_bar$labels,
    test_vignette_data$test_plot_core_dist_bar$labels
  )

  # Test line plot
  plot_core_dist_line <- plot_core_distribution(
    core_result = bcse_rare_core,
    core_set = "increase",
    group_var = "Crop",
    plot_type = "line"
  )

  expect_s3_class(plot_core_dist_line, "ggplot")
  expect_identical(
    plot_core_dist_line$labels,
    test_vignette_data$test_plot_core_dist_line$labels
  )

  # Test heatmap (with recoded Crop variable)
  bcse_rare_core$metadata <- bcse_rare_core$metadata |>
    mutate(
      Crop = recode(
        Crop,
        "Corn" = "Corn",
        "Sorghum" = "Sorghum",
        "Continuous Sorghum" = "Sorghum + cover crop",
        "Miscanthus" = "Miscanthus",
        "New Switchgrass" = "Establishing switchgrass",
        "Switchgrass" = "Mature switchgrass",
        "Early Succession" = "Successional vegetation",
        "Native Grasses" = "Native grass mix",
        "Prairie" = "Reconstructed prairie",
        "Poplar" = "Poplar"
      ),
      Crop = factor(
        Crop,
        levels = c(
          "Corn",
          "Sorghum",
          "Sorghum + cover crop",
          "Miscanthus",
          "Establishing switchgrass",
          "Mature switchgrass",
          "Successional vegetation",
          "Native grass mix",
          "Reconstructed prairie",
          "Poplar"
        )
      )
    )

  plot_core_dist_heatmap <- plot_core_distribution(
    core_result = bcse_rare_core,
    core_set = "increase",
    group_var = "Crop",
    plot_type = "heatmap"
  ) +
    viridis::scale_fill_viridis(option = "plasma", name = "Occupancy")

  expect_s3_class(plot_core_dist_heatmap, "ggplot")
  expect_identical(
    plot_core_dist_heatmap$labels$fill,
    test_vignette_data$test_plot_core_dist_heatmap$labels$fill
  )
})


test_that("4: Vignette neutral model fitting is consistent", {
  skip_on_cran()
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("tidyverse")

  # Load reference data
  load(here("tests/test_sets/test_vignette_data.rda"))

  # Recreate necessary objects
  data("bcse", package = "BRCore")
  library(phyloseq)
  library(tidyverse)

  # Use single thread on CI to ensure reproducibility
  n_threads <- if (nzchar(Sys.getenv("CI"))) {
    1
  } else if (.Platform$OS.type == "windows") {
    1
  } else {
    2
  }

  bcse_rarefied_otutable <- multi_rarefy(
    physeq = bcse,
    depth_level = 1000,
    num_iter = 100,
    threads = n_threads,
    set_seed = 7642
  )

  bcse_rare <- update_otu_table(
    physeq = bcse,
    otu_rare = bcse_rarefied_otutable
  )

  bcse_rare_core <- identify_core(
    physeq_obj = bcse_rare,
    priority_var = "Crop",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 2134
  )

  # Fit neutral model
  bcse_rare_core_neutral_fit <- fit_neutral_model(
    otu_table = bcse_rare_core$otu_table,
    core_set = bcse_rare_core$increase_core,
    abundance_occupancy = bcse_rare_core$abundance_occupancy
  )

  # Test: Neutral model fit structure
  expect_s3_class(bcse_rare_core_neutral_fit, "fit_neutral_model")
  expect_named(
    bcse_rare_core_neutral_fit,
    c("goodness_of_fit", "model_prediction")
  )

  # Test: Goodness of fit parameters match
  expect_equal(
    bcse_rare_core_neutral_fit$goodness_of_fit$m,
    test_vignette_data$test_bcse_rare_core_neutral_fit$goodness_of_fit$m,
    tolerance = 1e-2
  )

  expect_equal(
    bcse_rare_core_neutral_fit$goodness_of_fit$Rsqr,
    test_vignette_data$test_bcse_rare_core_neutral_fit$goodness_of_fit$Rsqr,
    tolerance = 1e-2
  )

  # Test: Model prediction dimensions match
  expect_equal(
    nrow(bcse_rare_core_neutral_fit$model_prediction),
    nrow(test_vignette_data$test_bcse_rare_core_neutral_fit$model_prediction)
  )

  # Test: Core membership classification matches
  core_membership_test <- bcse_rare_core_neutral_fit$model_prediction |>
    filter(membership == "Core")

  core_membership_ref <- test_vignette_data$test_bcse_rare_core_neutral_fit$model_prediction |>
    filter(membership == "Core")

  expect_equal(nrow(core_membership_test), nrow(core_membership_ref))
  expect_setequal(core_membership_test$otu, core_membership_ref$otu)

  # Plot neutral model
  plot_bcse_neutral_fit <- plot_neutral_model(bcse_rare_core_neutral_fit)

  # Test: Neutral model plot structure
  expect_s3_class(plot_bcse_neutral_fit, "ggplot")
  expect_identical(
    plot_bcse_neutral_fit$labels,
    test_vignette_data$test_plot_bcse_neutral_fit$labels
  )
})
