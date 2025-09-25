test_that("GlobalPatterns pipeline -> fit_neutral_model basic structure", {
    testthat::skip_on_cran()
    testthat::skip_if_not_installed("phyloseq")
    testthat::skip_if_not(exists("parallel_rarefy"),  "parallel_rarefy() not found")
    testthat::skip_if_not(exists("do_phyloseq"),      "do_phyloseq() not found")
    testthat::skip_if_not(exists("identify_core"),    "identify_core() not found")
    testthat::skip_if_not(exists("fit_neutral_model"),"fit_neutral_model() not found")
    
    # --- data & light prep (kept tiny for speed) ---
    data(GlobalPatterns, package = "phyloseq")
    testthat::expect_s4_class(GlobalPatterns, "phyloseq")
    
    rarefied_data <- parallel_rarefy(
        physeq      = GlobalPatterns,
        depth_level = 100,
        num_iter    = 3,
        threads     = 1,
        set_seed    = 999
    )
    
    rarefied_physeq <- do_phyloseq(
        physeq   = GlobalPatterns,
        otu_rare = rarefied_data
    )
    testthat::expect_s4_class(rarefied_physeq, "phyloseq")
    
    GP_core <- identify_core(
        physeq_obj       = rarefied_physeq,
        priority_var     = "SampleType",
        increase_value   = 0.02,
        abundance_weight = 0,
        seed             = 2324
    )
    
    # sanity on core output
    testthat::expect_true(all(c("otu_table", "increase_core", "occupancy_abundance") %in% names(GP_core)))
    
    # --- neutral model fit (matches how you call it) ---
    res <- fit_neutral_model(
        otu_table           = t(GP_core$otu_table),
        core_set            = GP_core$increase_core,
        abundance_occupancy = GP_core$occupancy_abundance
    )
    
    # basic structure only
    testthat::expect_type(res, "list")
    testthat::expect_true(all(c("model_prediction", "goodness_of_fit") %in% names(res)))
    testthat::expect_s3_class(res$model_prediction, "data.frame")
    testthat::expect_s3_class(res$goodness_of_fit,  "data.frame")
    testthat::expect_gt(nrow(res$model_prediction), 0)
    
    # OPTIONAL: quick plot smoke test (remove if you want this file to be faster)
    if (exists("plot_neutral_model")) {
        p <- plot_neutral_model(res)
        testthat::expect_s3_class(p, "ggplot")
        testthat::expect_silent(ggplot2::ggplot_build(p))
    }
})
