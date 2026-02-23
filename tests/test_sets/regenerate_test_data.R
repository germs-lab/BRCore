# Script to regenerate test_vignette_data.rda
# Run this script after making changes to multi_rarefy() or other functions
# that affect the vignette workflow outputs.
#
# Usage:
#   cd BRCore
#   Rscript tests/test_sets/regenerate_test_data.R
#
# This will regenerate the reference data used by test-vignette-workflow.R

library(BRCore)
library(phyloseq)
library(tidyverse)
library(viridis)

cat("=== Regenerating test_vignette_data.rda ===\n\n")

# Step 1: Load bcse data
cat("Loading bcse data...\n")
data("bcse", package = "BRCore")

# Step 2: Calculate rarefaction metrics
cat("Calculating rarefaction metrics...\n")
test_bcse_metrics <- add_rarefaction_metrics(data = bcse)

# Step 3: Plot rarefaction metrics
cat("Plotting rarefaction metrics...\n")
test_rarefaction_plot <- plot_rarefaction_metrics(test_bcse_metrics)

# Step 4: Multiple rarefaction
# Use single thread for consistency with test expectations
cat("Running multi_rarefy (this may take a while)...\n")
test_bcse_rarefied_otutable <- multi_rarefy(
  physeq = bcse,
  depth_level = 1000,
  num_iter = 100,
  threads = 1, # Single thread for consistency with CI test configuration
  set_seed = 7642
)

# Step 5: Update OTU table
cat("Updating OTU table...\n")
test_bcse_rare <- update_otu_table(
  physeq = bcse,
  otu_rare = test_bcse_rarefied_otutable
)

# Step 6: Identify core microbiome
cat("Identifying core microbiome...\n")
test_bcse_rare_core <- identify_core(
  physeq_obj = test_bcse_rare,
  priority_var = "Crop",
  increase_value = 0.02,
  abundance_weight = 0,
  seed = 2134
)

# Step 7: Plot identified core
cat("Plotting identified core...\n")
test_bcse_identified_core <- plot_identified_core(
  bray_curtis_ranked = test_bcse_rare_core$bray_curtis_ranked,
  elbow = test_bcse_rare_core$elbow,
  lastCall = test_bcse_rare_core$bc_increase,
  increase_value = test_bcse_rare_core$increase_value
)

# Step 8: Plot abundance-occupancy (increase)
cat("Plotting abundance-occupancy (increase)...\n")
test_plot_abund_occ_increase <- plot_abundance_occupancy(
  core_result = test_bcse_rare_core,
  core_set = "increase"
)

# Step 9: Plot abundance-occupancy (elbow)
cat("Plotting abundance-occupancy (elbow)...\n")
test_plot_abund_occ_elbow <- plot_abundance_occupancy(
  core_result = test_bcse_rare_core,
  core_set = "elbow"
)

# Step 10: Core distribution plots
cat("Creating core distribution plots...\n")

test_plot_core_dist_bar <- plot_core_distribution(
  core_result = test_bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "bar"
)

test_plot_core_dist_line <- plot_core_distribution(
  core_result = test_bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "line"
)

# Recode Crop variable for heatmap
test_bcse_rare_core_heatmap <- test_bcse_rare_core
test_bcse_rare_core_heatmap$metadata <- test_bcse_rare_core_heatmap$metadata |>
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

test_plot_core_dist_heatmap <- plot_core_distribution(
  core_result = test_bcse_rare_core_heatmap,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "heatmap"
) +
  viridis::scale_fill_viridis(option = "plasma", name = "Occupancy")

# Step 11: Fit neutral model
cat("Fitting neutral model...\n")
test_bcse_rare_core_neutral_fit <- fit_neutral_model(
  otu_table = test_bcse_rare_core$otu_table,
  core_set = test_bcse_rare_core$increase_core,
  abundance_occupancy = test_bcse_rare_core$abundance_occupancy
)

# Step 12: Plot neutral model
cat("Plotting neutral model...\n")
test_plot_bcse_neutral_fit <- plot_neutral_model(
  test_bcse_rare_core_neutral_fit
)

# Assemble the test data list
cat("\nAssembling test data...\n")
test_vignette_data <- list(
  test_bcse_metrics = test_bcse_metrics,
  test_rarefaction_plot = test_rarefaction_plot,
  test_bcse_rarefied_otutable = test_bcse_rarefied_otutable,
  test_bcse_rare = test_bcse_rare,
  test_bcse_rare_core = test_bcse_rare_core,
  test_bcse_identified_core = test_bcse_identified_core,
  test_plot_abund_occ_increase = test_plot_abund_occ_increase,
  test_plot_abund_occ_elbow = test_plot_abund_occ_elbow,
  test_plot_core_dist_bar = test_plot_core_dist_bar,
  test_plot_core_dist_line = test_plot_core_dist_line,
  test_plot_core_dist_heatmap = test_plot_core_dist_heatmap,
  test_bcse_rare_core_neutral_fit = test_bcse_rare_core_neutral_fit,
  test_plot_bcse_neutral_fit = test_plot_bcse_neutral_fit
)

# Save the test data
output_file <- "tests/test_sets/test_vignette_data.rda"
cat(paste0("Saving to ", output_file, "...\n"))
save(test_vignette_data, file = output_file)

cat("\n=== Done! ===\n")
cat("Reference data has been regenerated.\n")
cat("Please commit the updated test_vignette_data.rda file.\n")

# Print session info for reproducibility documentation
cat("\n=== Session Info ===\n")
print(sessionInfo())
