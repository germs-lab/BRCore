# Plot Bray-Curtis increase over ranked OTU/ASVs

Visualize the cumulative normalized mean Bray-Curtis increase returned
by
[identify_core()](http://www.germslab.org/BRCore/reference/identify_core.md),
over ranked OTU/ASVs and shows cutoff points for elbow percent increase
methods.

## Usage

``` r
plot_identified_core(
  bray_curtis_ranked,
  elbow,
  lastCall,
  increase_value = 0.02
)
```

## Arguments

- bray_curtis_ranked:

  A tibble as returned by `identify_core()$bray_curtis_ranked`.

- elbow:

  The number of OTU/ASVs identified by the elbow method (Integer).

- lastCall:

  The number of OTU/ASVs identified by the last percent Bray-Curtis
  increase method (Integer).

- increase_value:

  The percent increase value in decimal (e.g. 0.02) used for the
  Bray-Curtis increase method.

## Value

A ggplot2 object.

## Details

The function converts `rank` to integers and zooms the x-axis to the
first `1.2 * lastCall` ranks. Label positions are computed dynamically
from the observed `proportionBC` range to avoid overlap.

## See also

[identify_core()](http://www.germslab.org/BRCore/reference/identify_core.md)

[`identify_core`](http://www.germslab.org/BRCore/reference/identify_core.md)

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
# Example with the package switchgrass dataset
data("switchgrass", package = "BRCore")

# Identify core taxa
res <- identify_core(
  physeq_obj = switchgrass,
  priority_var = "sampling_date",
  increase_value = 0.02,
  seed = 48821
)
#> Seed used: 48821
#> ✔ Input phyloseq object is valid!
#> ℹ otu_table() is rarefied at a depth of: 1000
#> ℹ No taxonomy found (or empty). Continuing without taxonomy.
#> ✔ Core prioritizing variable: sampling_date
#> ℹ Ranked by Index only
#> ℹ Ranking OTUs based on BC dissimilarity, starting at 2026-04-21 16:22:33.751216
#> ✔ Elbow method identified 3 core OTUs
#> ✔ % increase method identified 34 core OTUs
#> ✔ Analysis complete!

# Plot using the returned curve and cut indices; label from increase_value
plot_identified_core(
  bray_curtis_ranked = res$bray_curtis_ranked,
  elbow = res$elbow,
  lastCall = res$bc_increase,
  increase_value = res$increase_value
)

# }
```
