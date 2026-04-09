# Test Sets Directory

This directory contains `test_vignette_data.rda` produced on 2026-02-05 by BAR. 

## Purpose

This file serves as the reference dataset for validating the complete BRCore workflow implementation. It is used when executing `tests/testthat/test-vignette-workflow.R`, which acts as the litmus test for the package's core analysis pipeline.

## Contents

`test_vignette_data.rda` is a list that contains saved outputs from a complete run of the BRCore vignette workflow, including:

- Rarefied phyloseq objects
- Core microbiome identification results (elbow and % increase methods)
- Abundance-occupancy data
- Neutral model fits
- All intermediate data structures and plots

## Usage

The test suite in `test-vignette-workflow.R` compares fresh workflow executions against these reference data to ensure:

- Reproducibility of results (with fixed seeds)
- Consistency of core identification algorithms
- Stability of plotting functions
- Integrity of the complete analysis pipeline

## Regenerating Reference Data

After making changes to `multi_rarefy()` or other functions that affect the vignette workflow outputs, the reference data must be regenerated to match the new implementation.

### Using the regeneration script

```bash
cd BRCore
Rscript tests/test_sets/regenerate_test_data.R
```

### Manual regeneration

Alternatively, run the complete vignette workflow:

```r
# Run the complete vignette workflow
knitr::purl("vignettes/BRCore-vignette.Rmd.orig", output = "vignettes/BRCore-vignette.R")
source("vignettes/BRCore-vignette.R")

# Save all intermediate objects
test_vignette_data <- list(
  test_bcse_metrics = bcse_metrics,
  test_rarefaction_plot = rarefaction_plot,
  test_bcse_rarefied_otutable = bcse_rarefied_list,
  test_bcse_rare = bcse_rare,
  test_bcse_rare_core = bcse_rare_core,
  test_bcse_identified_core = bcse_identified_core,
  test_plot_abund_occ_increase = plot_abund_occ_increase,
  test_plot_abund_occ_elbow = plot_abund_occ_elbow,
  test_plot_core_dist_bar = plot_core_dist_bar,
  test_plot_core_dist_line = plot_core_dist_line,
  test_plot_core_dist_heatmap = plot_core_dist_heatmap,
  test_bcse_rare_core_neutral_fit = bcse_rare_core_neutral_fit,
  test_plot_bcse_neutral_fit = plot_bcse_neutral_fit
)

save(test_vignette_data,
  file = "tests/test_sets/test_vignette_data.rda")