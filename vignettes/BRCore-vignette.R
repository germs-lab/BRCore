## ----include = FALSE----------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6,
  fig.path = "images/"
)
options(cli.progress_show_after = Inf) # Disable progress bar for cleaner output in vignette


## .note {
##   background-color: #e7f3fe;
##   border-left: 6px solid #2196F3;
##   padding: 10px;
##   margin: 15px 0;
## }

## ----setup, echo=TRUE---------------------------------------------------------------------------------
invisible(
  lapply(
    c("BRCore", "phyloseq", "tidyverse", "viridis"),
    library,
    character.only = TRUE
  )
)


## ----load bcse, echo=TRUE-----------------------------------------------------------------------------
data("bcse", package = "BRCore")
str(bcse)


## ----calculate metrics, echo=TRUE---------------------------------------------------------------------
bcse_metrics <- add_rarefaction_metrics(data = bcse)
bcse_metrics


## ----fig1_plot_metrics, echo=TRUE, fig.width=7, fig.height=6------------------------------------------
rarefaction_plot <- plot_rarefaction_metrics(bcse_metrics)
print(rarefaction_plot)


## ----rarefy bcse, echo=TRUE---------------------------------------------------------------------------
bcse_rarefied_list <-
  multi_rarefy(
    physeq_obj = bcse,
    depth_level = 1000,
    num_iter = 10,
    .as_array = FALSE,
    set_seed = 7642
  )


## ----verify success rarefaction, echo=TRUE------------------------------------------------------------
class(bcse_rarefied_list)
str(bcse_rarefied_list, max.level = 1)

## ----fig2_rarefaction_variance, echo=TRUE, fig.width=7, fig.height=6----------------------------------
rarefaction_variance_plot <- plot_variance_propagation(
  physeq_obj = bcse,
  rarefied = bcse_rarefied_list,
  q = 0,
  group_var = "Crop",
  group_color = "Plot"
)
rarefaction_variance_plot


## ----replace otu table, echo=TRUE---------------------------------------------------------------------
bcse_updated_rare <- update_otu_table(
  physeq_obj = bcse,
  rarefied_otus = bcse_rarefied_list,
  iteration = 1 # Speficify which iteration to use for the updated OTU table
)
print(bcse_updated_rare)
sample_sums(bcse_updated_rare)


## ----identify core microbiome, echo=TRUE--------------------------------------------------------------
bcse_rare_core <- identify_core(
  physeq_obj = bcse,
  priority_var = "Crop",
  increase_value = 0.02,
  abundance_weight = 0,
  depth_level = 1000,
  num_iter = 10,
  seed = 2134
)


## ----check the identified core, echo=TRUE-------------------------------------------------------------
str(bcse_rare_core)


## ----plot identified core bcse, echo=TRUE-------------------------------------------------------------
bcse_identified_core <- plot_identified_core(
  bray_curtis_ranked = bcse_rare_core$bray_curtis_ranked,
  elbow = bcse_rare_core$elbow,
  lastCall = bcse_rare_core$bc_increase,
  increase_value = bcse_rare_core$increase_value
)


## ----fig3_identified_core, echo=TRUE, fig.width=7, fig.height=6---------------------------------------
print(bcse_identified_core)


## ----plot abundance occupany and increase core set, echo=TRUE-----------------------------------------
plot_abund_occ_increase <- plot_abundance_occupancy(
  core_result = bcse_rare_core,
  core_set = "increase"
)


## ----fig4_plot_increase, echo=TRUE, fig.width=7, fig.height=6-----------------------------------------
print(plot_abund_occ_increase)


## ----plot abundance occupany and elbow core set, echo=TRUE--------------------------------------------
plot_abund_occ_elbow <- plot_abundance_occupancy(
  core_result = bcse_rare_core,
  core_set = "elbow"
)


## ----fig5_plot_elbow, echo=TRUE, fig.width=7, fig.height=6--------------------------------------------
plot_abund_occ_elbow +
  scale_fill_manual(values = c("darkgreen", "grey"))


## ----plot_type bar, echo=TRUE-------------------------------------------------------------------------
plot_core_dist_bar <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "bar"
)


## ----fig6_plot_type_bar, echo=TRUE--------------------------------------------------------------------
print(plot_core_dist_bar)


## ----plot_type line, echo=TRUE, fig.width=7, fig.height=6---------------------------------------------
plot_core_dist_line <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "line"
)


## ----fig7_plot_type_line, echo=TRUE-------------------------------------------------------------------
print(plot_core_dist_line)


## ----recode Crop for ease of plotting, echo=TRUE------------------------------------------------------
bcse_rare_core$metadata <- bcse_rare_core$metadata %>%
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


## ----plot_type heatmap, echo=TRUE, fig.width=7, fig.height=6------------------------------------------
plot_core_dist_heatmap <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "heatmap"
) +
  viridis::scale_fill_viridis(option = "plasma", name = "Occupancy")


## ----fig8_plot_type_heatmap, echo=TRUE----------------------------------------------------------------
print(plot_core_dist_heatmap)


## ----fit neutral model, echo=TRUE---------------------------------------------------------------------
bcse_rare_core_neutral_fit <- fit_neutral_model(
  otu_table = bcse_rare_core$otu_table,
  core_set = bcse_rare_core$increase_core,
  abundance_occupancy = bcse_rare_core$abundance_occupancy
)


## ----neutral fit result, echo=TRUE--------------------------------------------------------------------
str(bcse_rare_core_neutral_fit)


## ----plot neutral model, echo=TRUE--------------------------------------------------------------------
plot_bcse_neutral_fit <- plot_neutral_model(bcse_rare_core_neutral_fit)


## ----fig9_plot_neutral_fit, echo=TRUE, fig.width=7, fig.height=6--------------------------------------
print(plot_bcse_neutral_fit)


# Supplemental Information 

bcse_rare_core_iter1 <- identify_core(
    physeq_obj = bcse_updated_rare,
    priority_var = "Crop",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 2135
)


bcse_rare_core_iter1_w04 <- identify_core(
    physeq_obj = bcse_updated_rare,
    priority_var = "Crop",
    increase_value = 0.02,
    abundance_weight = 0.5,
    seed = 2135
)


## ----session info, echo=TRUE--------------------------------------------------------------------------
sessionInfo()
