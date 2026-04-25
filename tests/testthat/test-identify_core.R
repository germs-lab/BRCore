test_that("identify_core returns correctly formatted outputs on switchgrass", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  skip_if_not(is.function(identify_core), "identify_core() not found")

  # load test data shipped with the package
  data(
    "switchgrass",
    package = "BRCore",
    envir = environment()
  )
  skip_if_not(exists("switchgrass"), "switchgrass dataset not found")

  # sanity: required grouping column
  sm <- as(phyloseq::sample_data(switchgrass), "data.frame")
  skip_if_not(
    "sampling_date" %in% names(sm),
    "'sampling_date' column not found in sample_data(switchgrass)"
  )

  # run once with a fixed seed for reproducibility
  res <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 4
  )

  # ---- presence / type checks ----
  expect_true(is.list(res))
  expect_true(all(
    c(
      "bray_curtis_ranked",
      "otu_ranked",
      "abundance_occupancy",
      "elbow",
      "bc_increase"
    ) %in%
      names(res)
  ))

  # bray_curtis_ranked
  expect_true(is.data.frame(res$bray_curtis_ranked))
  expect_true(all(
    c("rank", "MeanBC", "proportionBC", "IncreaseBC") %in%
      names(res$bray_curtis_ranked)
  ))
  expect_true(is.numeric(res$bray_curtis_ranked$MeanBC))
  expect_true(is.numeric(res$bray_curtis_ranked$proportionBC))
  expect_true(is.numeric(res$bray_curtis_ranked$IncreaseBC))
  expect_equal(max(res$bray_curtis_ranked$proportionBC), 1)
  expect_true(is.na(res$bray_curtis_ranked$IncreaseBC[1]))

  # otu_ranked
  expect_true(is.data.frame(res$otu_ranked))
  expect_true(all(c("otu", "rank") %in% names(res$otu_ranked)))
  expect_true(is.numeric(res$otu_ranked$rank))
  expect_gt(nrow(res$otu_ranked), 0)

  # abundance occupancy
  expect_true(is.data.frame(res$abundance_occupancy))
  expect_true(all(
    c("otu", "otu_occ", "otu_rel") %in% names(res$abundance_occupancy)
  ))
  expect_true(is.numeric(res$abundance_occupancy$otu_occ))
  expect_true(is.numeric(res$abundance_occupancy$otu_rel))
  expect_true(all(
    res$abundance_occupancy$otu_occ >= 0 &
      res$abundance_occupancy$otu_occ <= 1
  ))
  expect_true(all(
    is.finite(res$abundance_occupancy$otu_rel) &
      res$abundance_occupancy$otu_rel >= 0
  ))
  # same OTU universe (ignoring order)
  expect_true(setequal(res$abundance_occupancy$otu, res$otu_ranked$otu))

  # elbow & bc_increase are length-1 integers/numerics within bounds
  expect_length(as.integer(res$elbow), 1)
  expect_length(as.integer(res$bc_increase), 1)
  expect_gte(as.integer(res$elbow), 1)
  expect_gte(as.integer(res$bc_increase), 1)
  expect_lte(as.integer(res$elbow), nrow(res$otu_ranked))
  expect_lte(as.integer(res$bc_increase), nrow(res$otu_ranked))

  # optional: number of BC ranks should match number of taxa considered
  expect_equal(nrow(res$bray_curtis_ranked), nrow(res$otu_ranked))

  # ---- reproducibility with same seed ----
  res2 <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 4
  )
  expect_equal(
    res$bray_curtis_ranked$MeanBC,
    res2$bray_curtis_ranked$MeanBC,
    tolerance = 1e-12
  )
  expect_equal(res$elbow, res2$elbow)
  expect_equal(res$bc_increase, res2$bc_increase)
})

test_that("max_otus parameter correctly limits the number of OTUs analyzed", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")

  data(
    "switchgrass",
    package = "BRCore",
    envir = environment()
  )
  skip_if_not(exists("switchgrass"), "switchgrass dataset not found")

  # Get total number of OTUs
  total_otus <- phyloseq::ntaxa(switchgrass)

  # Test 1: max_otus smaller than total OTUs should truncate
  max_test <- 50
  res_limited <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    max_otus = max_test,
    seed = 123
  )

  expect_equal(nrow(res_limited$otu_ranked), max_test)
  expect_equal(nrow(res_limited$bray_curtis_ranked), max_test)
  expect_lte(res_limited$elbow, max_test)
  expect_lte(res_limited$bc_increase, max_test)

  # Test 2: max_otus larger than total OTUs should use all
  res_all <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    max_otus = total_otus + 100,
    seed = 123
  )

  expect_equal(nrow(res_all$otu_ranked), total_otus)
  expect_equal(nrow(res_all$bray_curtis_ranked), total_otus)

  # Test 3: NULL max_otus should use all OTUs
  res_null <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    max_otus = NULL,
    seed = 123
  )

  expect_equal(nrow(res_null$otu_ranked), total_otus)
  expect_equal(nrow(res_null$bray_curtis_ranked), total_otus)

  # Test 4: Results with NULL and large max_otus should be identical
  expect_equal(res_all$otu_ranked$otu, res_null$otu_ranked$otu)
  expect_equal(
    res_all$bray_curtis_ranked$MeanBC,
    res_null$bray_curtis_ranked$MeanBC
  )
  expect_equal(res_all$elbow, res_null$elbow)
  expect_equal(res_all$bc_increase, res_null$bc_increase)
})

test_that("max_otus parameter validates inputs correctly", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")

  data(
    "switchgrass",
    package = "BRCore",
    envir = environment()
  )
  skip_if_not(exists("switchgrass"), "switchgrass dataset not found")

  # Test invalid inputs
  expect_error(
    identify_core(
      physeq_obj = switchgrass,
      priority_var = "sampling_date",
      max_otus = -5,
      seed = 123
    ),
    "positive integer"
  )

  expect_error(
    identify_core(
      physeq_obj = switchgrass,
      priority_var = "sampling_date",
      max_otus = 0,
      seed = 123
    ),
    "positive integer"
  )

  expect_error(
    identify_core(
      physeq_obj = switchgrass,
      priority_var = "sampling_date",
      max_otus = "fifty",
      seed = 123
    ),
    "positive integer"
  )
})

test_that("max_otus preserves top-ranked OTUs based on ranking method", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")

  data(
    "switchgrass",
    package = "BRCore",
    envir = environment()
  )
  skip_if_not(exists("switchgrass"), "switchgrass dataset not found")

  # Run with all OTUs
  res_full <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    seed = 456
  )

  # Run with limited OTUs
  max_test <- 30
  res_limited <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    max_otus = max_test,
    seed = 456
  )

  # The limited set should contain the top N OTUs from the full ranking
  top_otus_full <- head(res_full$otu_ranked$otu, max_test)

  expect_equal(res_limited$otu_ranked$otu, top_otus_full)
  expect_true(all(res_limited$otu_ranked$otu %in% res_full$otu_ranked$otu))
})

test_that("max_otus with abundance_weight preserves correct ranking", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")

  data(
    "switchgrass",
    package = "BRCore",
    envir = environment()
  )
  skip_if_not(exists("switchgrass"), "switchgrass dataset not found")

  # Test with abundance weighting
  res_full_weighted <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    abundance_weight = 0.3,
    seed = 789
  )

  max_test <- 40
  res_limited_weighted <- identify_core(
    physeq_obj = switchgrass,
    priority_var = "sampling_date",
    abundance_weight = 0.3,
    max_otus = max_test,
    seed = 789
  )

  # Check that the limited set matches the top N from full weighted ranking
  top_weighted <- head(res_full_weighted$otu_ranked$otu, max_test)
  expect_equal(res_limited_weighted$otu_ranked$otu, top_weighted)
  expect_equal(nrow(res_limited_weighted$otu_ranked), max_test)
})

# Test identify_core internal functions if accesible
