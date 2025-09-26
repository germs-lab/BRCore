test_that("GlobalPatterns pipeline -> fit_neutral_model basic structure", {
    testthat::skip_on_cran()
    testthat::skip_if_not_installed("phyloseq")
    testthat::skip_if_not(exists("parallel_rarefy"),   "parallel_rarefy() not found")
    testthat::skip_if_not(exists("do_phyloseq"),       "do_phyloseq() not found")
    testthat::skip_if_not(exists("identify_core"),     "identify_core() not found")
    testthat::skip_if_not(exists("fit_neutral_model"), "fit_neutral_model() not found")
    
    # make partial-arg bugs loud
    old <- getOption("warnPartialMatchArgs"); options(warnPartialMatchArgs = TRUE)
    on.exit(options(warnPartialMatchArgs = old), add = TRUE)
    
    data(GlobalPatterns, package = "phyloseq")
    expect_s4_class(GlobalPatterns, "phyloseq")
    
    # 1) prune to depth first
    target_depth <- 1000  # use 500 if you prefer, but keep prune consistent with it
    GP <- phyloseq::prune_samples(phyloseq::sample_sums(GlobalPatterns) >= target_depth, GlobalPatterns)
    
    # 2) rarefy
    rarefied_otu <- parallel_rarefy(
        physeq      = GP,
        depth_level = target_depth,
        num_iter    = 3,
        threads     = 1,
        set_seed    = 991
    )
    
    # 3) rebuild phyloseq with the EXACT arg name your function expects
    GlobalPatterns_rare <- do_phyloseq(
        physeq        = GP,
        otu_rare = rarefied_otu   # <-- change this name if your function uses a different one
    )
    expect_s4_class(GlobalPatterns_rare, "phyloseq")
    
    # 4) assert even depth before identify_core()
    depths <- phyloseq::sample_sums(GlobalPatterns_rare)
    expect_true(length(unique(depths)) == 1)
    
    GlobalPatterns_core <- identify_core(
        physeq_obj       = GlobalPatterns_rare,
        priority_var     = "SampleType",
        increase_value   = 0.02,
        abundance_weight = 0,
        seed             = 123
    )
    expect_true(all(c("otu_table", "increase_core", "occupancy_abundance") %in% names(GlobalPatterns_core)))
    
    GlobalPatterns_fit <- fit_neutral_model(
        otu_table           = t(GlobalPatterns_core$otu_table),
        core_set            = GlobalPatterns_core$increase_core,
        abundance_occupancy = GlobalPatterns_core$occupancy_abundance
    )
    
    expect_type(GlobalPatterns_fit, "list")
    expect_true(all(c("model_prediction", "goodness_of_fit") %in% names(GlobalPatterns_fit)))
    expect_s3_class(GlobalPatterns_fit$model_prediction, "data.frame")
    expect_s3_class(GlobalPatterns_fit$goodness_of_fit,  "data.frame")
    expect_gt(nrow(GlobalPatterns_fit$model_prediction), 0)
    
    if (exists("plot_neutral_model")) {
        p <- plot_neutral_model(GlobalPatterns_fit)
        expect_s3_class(p, "ggplot")
        expect_silent(ggplot2::ggplot_build(p))
    }
})
