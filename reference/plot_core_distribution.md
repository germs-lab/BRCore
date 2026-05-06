# Plot Core Taxa Occupancy Across Metadata Groups

Creates a plot showing core taxa (i.e. OTUs/ASVs) occupancy patterns
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
  [`identify_core`](https://www.germslab.org/BRCore/reference/identify_core.md),
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

## See also

[`identify_core()`](https://www.germslab.org/BRCore/reference/identify_core.md),
[`plot_abundance_occupancy()`](https://www.germslab.org/BRCore/reference/plot_abundance_occupancy.md),
and
[`plot_identified_core()`](https://www.germslab.org/BRCore/reference/plot_identified_core.md)

## Examples

``` r
library(BRCore)
data("switchgrass_core", package = "BRCore")

p <- plot_core_distribution(
  core_result = switchgrass_core,
  core_set = "increase",
  group_var = "sampling_date",
  plot_type = "bar"
)
print(p)

```
