# Plot Core Taxa Occupancy Across Metadata Groups

Creates a bar plot showing core taxa (i.e. OTUs/ASVs) occupancy patterns
across a grouping variable.

## Usage

``` r
plot_core_distribution(
  core_result,
  core_set = "elbow",
  group_var = "Crop",
  plot_type = c("bar", "line", "heatmap")
)
```

## Arguments

- core_result:

  A list object returned by
  [`identify_core`](http://www.germslab.org/BRCore/reference/identify_core.md),
  containing at minimum:

  - `otu_table`: A data frame with ASV/OTUs as rows and samples as
    columns.

  - `metadata`: A data frame with samples as rows and grouping variables
    as columns.

  - `otu_ranked`: A data frame with ranked taxa containing:

    - `otu`: A column with taxa names.

    - `rank`: A column with the rank for each taxon.

  - `elbow_core`: Character vector of OTU IDs identified as core using
    the elbow method.

  - `increase_core`: Character vector of OTU IDs identified as core
    using the increase method.

- core_set:

  Which core set to plot: "elbow" (default) or "increase".

- group_var:

  Metadata column for bar coloring. Default: "sampling_date".

- plot_type:

  Allows selection of 3 different plot types: `bar`, `line`, or
  `heatmap`.

## Value

A ggplot2 object that can be further customized.

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
# Generate an object from identify_core and then plot
data("switchgrass", package = "BRCore")

switchgrass_core <- identify_core(
  physeq_obj = switchgrass,
  priority_var = "sampling_date",
  increase_value = 0.02,
  abundance_weight = 0,
  seed = 1234
)
#> Seed used: 1234
#> ✔ Input phyloseq object is valid!
#> ℹ No `rarefied_list` provided. `physeq_obj` is already rarefied; wrapping as a single iteration.
#> ℹ No taxonomy found (or empty). Continuing without taxonomy.
#> ✔ Core prioritizing variable: sampling_date
#> ℹ Ranked by Rank only
#> ℹ Ranking OTUs based on BC dissimilarity, starting at 2026-04-25 19:01:40.870299
#> ■■■■■                             15% | ETA:  6s
#> ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ■■■■■■■■■■■■■■■■■■■■■             67% | ETA:  3s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% | ETA:  1s
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
#> ✔ Elbow method identified 3 core OTUs
#> ✔ % increase method identified 34 core OTUs
#> ✔ Analysis complete!

plot_core_distribution(
  core_result = switchgrass_core,
  core_set = "increase",
  group_var = "sampling_date",
  plot_type = "bar"
)

# }
```
