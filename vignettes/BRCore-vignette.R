## ----chunk setup, include = FALSE--------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
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
## 
## h1, .h1 {
##     margin-top: 84px;
##     margin-bottom: 42px;
## }
## h2, .h2, h3, .h3 {
##     margin-top: 42px;
##     margin-bottom: 21px;
## }
## 
## p.caption {
##     font-size: 1em;
##     font-style: italic;
##     color: grey;
##     margin-right: 10%;
##     margin-left: 10%;
##     text-align: justify;
## }

## ----load libraries, echo=TRUE-----------------------------------------------------------
invisible(
  lapply(
    c("BRCore", "phyloseq", "tidyverse", "viridis"),
    library,
    character.only = TRUE
  )
)


## ----load bcse, echo=TRUE----------------------------------------------------------------
data("bcse", package = "BRCore")
str(bcse)


## ----calculate metrics, echo=TRUE--------------------------------------------------------
bcse_metrics <- add_rarefaction_metrics(data = bcse)
bcse_metrics


## ----fig1_plot_metrics, echo=TRUE, fig.cap="Figure 1: Rarefaction metrics. a-b, histograms of samples. c, Good's coverage per total number of sequences in a sample. d-e, log10 of sequences in a sample. f,  samples ranked by sequence reads in a sample."----
rarefaction_plot <- plot_rarefaction_metrics(bcse_metrics)
print(rarefaction_plot)


## ----rarefy bcse, echo=TRUE--------------------------------------------------------------
bcse_rarefied_list <-
  multi_rarefy(
    physeq_obj = bcse,
    depth_level = 1000,
    num_iter = 10,
    .as_array = FALSE,
    set_seed = 7642
  )


## ----verify success rarefaction, echo=TRUE-----------------------------------------------
class(bcse_rarefied_list)
str(bcse_rarefied_list[[1]], list.len = 10) # Dimensions of iteration #1


## ----fig2_rarefaction_variance, echo=TRUE, fig.cap="Figure 2: Alpha diversity variance propagation (variance between rarefaction iteration) across 10 iterations of rarefaction. Points represent samples. The  x-axis represent the grouping variable (Crop). The y-axis represent the alpha diversity metric (q = 0, i.e. richness) calculated on the samples. The color of the points represent the 'Plot' variable, which is a nested variable within 'Crop'."----
rarefaction_variance_plot <- plot_variance_propagation(
  physeq_obj = bcse,
  rarefied = bcse_rarefied_list,
  q = 0,
  group_var = "Crop",
  group_color = "Plot"
)
print(rarefaction_variance_plot)


## ----update_otu_table, echo=TRUE---------------------------------------------------------
bcse_updated_rare <- update_otu_table(
  physeq_obj = bcse,
  rarefied_otus = bcse_rarefied_list,
  iteration = 1 # Speficify which iteration to use for the updated OTU table
)


## ----identify_core, echo=TRUE------------------------------------------------------------
bcse_rare_core <- identify_core(
  physeq_obj = bcse,
  priority_var = "Crop",
  increase_value = 0.02,
  abundance_weight = 0,
  depth_level = 1000,
  num_iter = 10,
  seed = 2134
)


## ----check the identified core, echo=TRUE------------------------------------------------
str(bcse_rare_core)


## ----fig3_identified_core, echo=TRUE, fig.cap="Figure 3: Percent Bray-Curtis similarity between samples per ranked ASV/OTUs. Number of core ASV/OTUs identified by Elbow and Last 2% increase in Bray-Curtis similarity are shown."----
bcse_identified_core <- plot_identified_core(
  bray_curtis_ranked = bcse_rare_core$bray_curtis_ranked,
  elbow = bcse_rare_core$elbow,
  lastCall = bcse_rare_core$bc_increase,
  increase_value = bcse_rare_core$increase_value
)

print(bcse_identified_core)


## ----plot abundance occupany and increase core set, echo=TRUE----------------------------


## ----fig4_plot_increase, echo=TRUE, fig.cap="Figure 4: Abundance-occupancy distribution for the 'bcse' dataset. The core ASV/OTUs identified by the last 2% increase method are highlighted in red."----
plot_abund_occ_increase <- plot_abundance_occupancy(
  core_result = bcse_rare_core,
  core_set = "increase"
)
print(plot_abund_occ_increase)


## ----fig5_plot_elbow, echo=TRUE, fig.cap="Figure 5: Abundance-occupancy distribution for the 'bcse' dataset. The core ASV/OTUs identified by the elbow method are highlighted in dark green."----
plot_abund_occ_elbow <- plot_abundance_occupancy(
  core_result = bcse_rare_core,
  core_set = "elbow"
)

plot_abund_occ_elbow +
  scale_fill_manual(values = c("darkgreen", "grey"))


## ----fig6_plot_type_bar, echo=TRUE, fig.cap="Figure 6: Occupancy of core ASV/OTUs across the 'Crop' variable. Each bar represents the average occupancy of core ASV/OTUs in samples belonging to each level of the 'Crop' variable. A bar plot won't work here, way too many variable levels."----
plot_core_dist_bar <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "bar"
)

print(plot_core_dist_bar)


## ----plot_type line----------------------------------------------------------------------



## ----fig7_plot_type_line, echo=TRUE, fig.width=7, fig.height=10, fig.cap="Figure 7: Occupancy of core ASV/OTUs across the 'Crop' variable. Each point represents the average occupancy of core ASV/OTUs in samples belonging to each level of the 'Crop' variable. A line plot is better than a bar plot but still not ideal for this many variable levels."----
plot_core_dist_line <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "line"
)
print(plot_core_dist_line)


## ----reorder variable levels, echo=TRUE--------------------------------------------------
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


## ----plot_type heatmap-------------------------------------------------------------------



## ----fig8_plot_type_heatmap, echo=TRUE, fig.cap="Figure 8: Occupancy of core ASV/OTUs across the 'Crop' variable. Each cell represents the average occupancy of core ASV/OTUs in samples belonging to each level of the 'Crop' variable. A heatmap is more compact and shows well enough the average occupancy across samples in each variable level."----
plot_core_dist_heatmap <- plot_core_distribution(
  core_result = bcse_rare_core,
  core_set = "increase",
  group_var = "Crop",
  plot_type = "heatmap"
) +
  viridis::scale_fill_viridis(option = "plasma", name = "Occupancy")

print(plot_core_dist_heatmap)


## ----fit neutral model, echo=TRUE--------------------------------------------------------
bcse_rare_core_neutral_fit <- fit_neutral_model(
  otu_table = bcse_rare_core$otu_table,
  core_set = bcse_rare_core$increase_core,
  abundance_occupancy = bcse_rare_core$abundance_occupancy
)


## ----neutral fit result, echo=TRUE-------------------------------------------------------
str(bcse_rare_core_neutral_fit)


## ----fig9_plot_neutral_fit, echo=TRUE, fig.cap="Figure 9: Neutral model fit illustrates the neutral model of abundance-occupancy distributions for the 'bcse' dataset. The $R^2$ value represents a standard coefficient of determination, calculated as $R^2 = 1 - \\frac{SS_{err}}{SS_{total}}$, providing a measure of goodness of fit. The immigration parameter ($m$) estimates the probability that an individual in a local sample originates from the metacommunity rather than from local dispersal. We also observe that the neutral model provides a poor fit; the $R^2$ is notably low, and the majority of taxa (including both core and non-core ASV/OTUs) fall above the model's prediction line. Furthermore, the estimated immigration parameter ($m = 0.63$) suggests that the local community is dominated by strong dispersal and mixing from the metacommunity. "----
plot_bcse_neutral_fit <- plot_neutral_model(bcse_rare_core_neutral_fit)

print(plot_bcse_neutral_fit)


## ----supplemental_info, echo=TRUE--------------------------------------------------------
bcse_rare_core_iter1 <- identify_core(
  physeq_obj = bcse_updated_rare,
  priority_var = "Crop",
  increase_value = 0.02,
  abundance_weight = 0,
  seed = 2135
)

