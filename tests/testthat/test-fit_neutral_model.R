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
    
    rarefied_otu <- parallel_rarefy(
        physeq      = GlobalPatterns,
        depth_level = 100,
        num_iter    = 3,
        threads     = 1,
        set_seed    = 999
    )
    
    GlobalPatterns_rare <- do_phyloseq(
        physeq   = GlobalPatterns,
        otu_rare = rarefied_otu
    )
    testthat::expect_s4_class(GlobalPatterns_rare, "phyloseq")
    
    GlobalPatterns_core <- identify_core(
        physeq_obj       = GlobalPatterns_rare,
        priority_var     = "SampleType",
        increase_value   = 0.02,
        abundance_weight = 0,
        seed             = 123
    )
    
    # sanity on core output
    testthat::expect_true(all(c("otu_table", "increase_core", "occupancy_abundance") %in% names(GlobalPatterns_core)))
    
    # --- neutral model fit (matches how you call it) ---
    GlobalPatterns_fit <- fit_neutral_model(
        otu_table           = t(GlobalPatterns_core$otu_table),
        core_set            = GlobalPatterns_core$increase_core,
        abundance_occupancy = GlobalPatterns_core$occupancy_abundance
    )
    
    # basic structure only
    testthat::expect_type(GlobalPatterns_fit, "list")
    testthat::expect_true(all(c("model_prediction", "goodness_of_fit") %in% names(GlobalPatterns_fit)))
    testthat::expect_s3_class(GlobalPatterns_fit$model_prediction, "data.frame")
    testthat::expect_s3_class(GlobalPatterns_fit$goodness_of_fit,  "data.frame")
    testthat::expect_gt(nrow(GlobalPatterns_fit$model_prediction), 0)
    
    # OPTIONAL: quick plot smoke test (remove if you want this file to be faster)
    if (exists("plot_neutral_model")) {
        p <- plot_neutral_model(GlobalPatterns_fit)
        testthat::expect_s3_class(p, "ggplot")
        testthat::expect_silent(ggplot2::ggplot_build(p))
    }
})
