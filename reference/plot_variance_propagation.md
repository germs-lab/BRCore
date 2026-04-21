# Variance propagation diagnostic for rarefaction

This function evaluate the variance generated during multiple
rarefaction by plotting comparing raw vs. rarefied alpha diversity
metrics calculated at each iterations. It is possible to plot observed
richness (q=0), Shannon diversity (q=1), or Simpson diversity (q=2) by
setting the `q` parameter to "richness" or `q` = 0, "shannon" or `q` =
1, or "shannon" or `q` = 2. The plot is faceted by method (raw vs
rarefied) and colored by a specified grouping variable from the sample
data.

## Usage

``` r
plot_variance_propagation(physeq_obj, rarefied, q = 0, group_var, group_color)
```

## Arguments

- physeq_obj:

  Raw phyloseq object

- rarefied:

  Output from multi_rarefy(). Either a list of dataframes or and array.

- q:

  Hill number order (q = 0 for richness, q = 1 for Shannon, q = 2 for
  Simpson)

- group_var:

  A grouping variable to use gor grouping as in the sample_data()

- group_color:

  A color variable to use present in the sample_data()

## Value

ggplot object comparing raw vs rarefied diversity distributions across
iterations.

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
# Example using the bcse dataset, comparing hill q=1 between Poplar and Switchgrass plots
bcse_filt <- bcse |>
subset_samples(Crop %in% c("Poplar", "Switchgrass"))
bcse_rarefied_otutable_filt <-
 multi_rarefy(
       physeq_obj = bcse_filt,
       depth_level = 1000,
       num_iter = 100,
       .as_array = FALSE,
       set_seed = 7643
   )
#> Error in multi_rarefy(physeq_obj = bcse_filt, depth_level = 1000, num_iter = 100,     .as_array = FALSE, set_seed = 7643): unused arguments (physeq_obj = bcse_filt, .as_array = FALSE)

plot_variance_propagation(
   physeq_obj   = bcse_filt,
   rarefied = bcse_rarefied_otutable_filt,
   q        = 1,
   group_var = "Crop",
   group_color = "Plot"
)
#> Error in plot_variance_propagation(physeq_obj = bcse_filt, rarefied = bcse_rarefied_otutable_filt,     q = 1, group_var = "Crop", group_color = "Plot"): could not find function "plot_variance_propagation"
# }
```
