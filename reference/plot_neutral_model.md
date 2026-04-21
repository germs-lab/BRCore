# Plot a fitted Neutral Model to Microbial Community Data

This function plots ASV/OTUs log abundances into a fitted neutral model
of microbial abundance-occupancy distribution. ASV/OTUs that are *Core*,
*As predicted* (i.e. Neutral), *Below* and *Above* model predictions are
drawn with distinct point colors, see details.

## Usage

``` r
plot_neutral_model(fit_result)
```

## Arguments

- fit_result:

  A list-like object returned by
  [fit_neutral_model()](http://www.germslab.org/BRCore/reference/fit_neutral_model.md).

## Value

A [ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object
with:

- *x*-axis: `log10(mean abundance)` (`log10(otu_rel)`);

- *y*-axis: `Occupancy` (`otu_occ`);

- Points colored by membership/fit class as described above;

- Neutral-model curve (solid) with 95\\

- An inset white label reporting *R*\\^2\\ and *m* taken directly from
  `fit_result$goodness_of_fit`.

The function does not recompute statistics; it only visualizes the
supplied predictions and metrics.

## Details

Points are split into four groups for display:

- **Not core (as predicted)** – background points;

- **Core (as predicted)** – core taxa whose occupancy matches the model;

- **Core (above prediction)** – core taxa above the 95\\

- **Core (below prediction)** – core taxa below the 95\\

## See also

[`fit_neutral_model`](http://www.germslab.org/BRCore/reference/fit_neutral_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have run a neutral model fit on the switchgrass data.frames
is.data.frame(switchgrass_core_fit$model_prediction)
is.data.frame(switchgrass_core_fit$goodness_of_fit)

p <- plot_neutral_model(switchgrass_core_fit)
print(p)
} # }
```
