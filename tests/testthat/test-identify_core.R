test_that("identify_core returns correctly formatted outputs on switchgrass", {
    skip_if_not_installed("phyloseq")
    skip_if_not_installed("vegan")
    skip_if_not(is.function(identify_core), "identify_core() not found")
    
    # load test data shipped with the package
    suppressWarnings(data("switchgrass", envir = environment()))
    skip_if_not(exists("switchgrass"), "switchgrass dataset not found")
    
    # sanity: required grouping column
    sm <- as(phyloseq::sample_data(switchgrass), "data.frame")
    skip_if_not("sampling_date" %in% names(sm),
                "'sampling_date' column not found in sample_data(switchgrass)")
    
    # run once with a fixed seed for reproducibility
    res <- identify_core(
        physeq_obj = switchgrass,
        priority_var = "sampling_date",
        increase_value = 0.02,
        seed          = 4
    )
    
    # ---- presence / type checks ----
    expect_true(is.list(res))
    expect_true(all(c("bray_curtis_ranked",
                      "otu_ranked",
                      "occupancy_abundance",
                      "elbow",
                      "bc_increase") %in% names(res)))
    
    # bray_curtis_ranked
    expect_true(is.data.frame(res$bray_curtis_ranked))
    expect_true(all(c("rank", "MeanBC", "proportionBC", "IncreaseBC") %in% names(res$bray_curtis_ranked)))
    expect_true(is.numeric(res$bray_curtis_ranked$MeanBC))
    expect_true(is.numeric(res$bray_curtis_ranked$proportionBC))
    expect_true(is.numeric(res$bray_curtis_ranked$IncreaseBC))
    expect_equal(max(res$bray_curtis_ranked$proportionBC), 1)
    expect_equal(res$bray_curtis_ranked$IncreaseBC[1], 0)
    
    # otu_ranked
    expect_true(is.data.frame(res$otu_ranked))
    expect_true(all(c("otu", "rank") %in% names(res$otu_ranked)))
    expect_true(is.numeric(res$otu_ranked$rank))
    expect_gt(nrow(res$otu_ranked), 0)
    
    # occupancy_abundance
    expect_true(is.data.frame(res$occupancy_abundance))
    expect_true(all(c("otu", "otu_occ", "otu_rel") %in% names(res$occupancy_abundance)))
    expect_true(is.numeric(res$occupancy_abundance$otu_occ))
    expect_true(is.numeric(res$occupancy_abundance$otu_rel))
    expect_true(all(res$occupancy_abundance$otu_occ >= 0 & res$occupancy_abundance$otu_occ <= 1))
    expect_true(all(is.finite(res$occupancy_abundance$otu_rel) & res$occupancy_abundance$otu_rel >= 0))
    # same OTU universe (ignoring order)
    expect_true(setequal(res$occupancy_abundance$otu, res$otu_ranked$otu))
    
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
        physeq_obj    = switchgrass,
        priority_var  = "sampling_date",
        increase_value = 0.02,
        abundance_weight = 0,
        seed = 4
    )
    expect_equal(res$bray_curtis_ranked$MeanBC, res2$bray_curtis_ranked$MeanBC, tolerance = 1e-12)
    expect_equal(res$elbow, res2$elbow)
    expect_equal(res$bc_increase, res2$bc_increase)
})
