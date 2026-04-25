# BRCore

## Overview

BRCore provides a unified framework for identification and ecological
interpretation of core microbiomes across time and space, enhancing
robustness and reproducibility in microbiome data analysis.

BRCore provides tools for:

**Rarefaction analysis**: Calculate pre-rarefaction metrics and perform
multiple rarefactions

- [`add_rarefaction_metrics()`](http://www.germslab.org/BRCore/reference/add_rarefaction_metrics.md)
- [`multi_rarefy()`](http://www.germslab.org/BRCore/reference/multi_rarefy.md)

**Core microbiome identification**: Identify core microbial taxa using
abundance-occupancy distributions

- [`identify_core()`](http://www.germslab.org/BRCore/reference/identify_core.md)

**Neutral model fitting**: Fit and visualize neutral community models

- [`fit_neutral_model()`](http://www.germslab.org/BRCore/reference/fit_neutral_model.md)
- [`plot_neutral_model()`](http://www.germslab.org/BRCore/reference/plot_neutral_model.md)

**Visualization**: Plot rarefaction diagnostics, abundance-occupancy
curves, and core distributions

- [`plot_rarefaction_metrics()`](http://www.germslab.org/BRCore/reference/plot_rarefaction_metrics.md)
- [`plot_abundance_occupancy()`](http://www.germslab.org/BRCore/reference/plot_abundance_occupancy.md)
- [`plot_core_distribution()`](http://www.germslab.org/BRCore/reference/plot_core_distribution.md)
- [`plot_identified_core()`](http://www.germslab.org/BRCore/reference/plot_identified_core.md)

## Installation

Install the latest *stable* version of BRCore from CRAN with:

``` r
install.packages("BRCore")
```

Install the *development* version of BRCore from GitHub with:

``` r
# install.packages("pak")
pak::pak("germs-lab/BRCore")
```

***Note:*** *`pak` handles dependencies automatically.*

## Quick Start

``` r
library(BRCore)
library(phyloseq)

# Load example data
data("bcse", package = "BRCore")

# Add rarefaction metrics
bcse_metrics <- add_rarefaction_metrics(data = bcse)

# Perform multiple rarefaction
bcse_rarefied_list <- multi_rarefy(
  physeq_obj = bcse,
  depth_level = 1000,
  num_iter = 3,
  set_seed = 7642
)

# Update phyloseq object with rarefied data
bcse_rare_single <- update_otu_table(physeq_obj = bcse, rarefied_otus = bcse_rarefied_list, iteration = 2) # Your preffered iteration can be used here

# Identify core microbiome

# With a single iteration of rarefaction
bcse_core <- identify_core(
  physeq_obj = bcse_rare_single,
  priority_var = "Crop",
  increase_value = 0.02,
  seed = 2134
)


# With multiple iterations of rarefaction
bcse_core_multi <- identify_core(
  physeq_obj = bcse, 
  rarefied_list = bcse_rarefied_list
  priority_var = "Crop",
  increase_value = 0.02,
  depth_level = 1000,
  num_iter = 10,
  seed = 2134
)

# Visualize abundance-occupancy distribution
plot_abundance_occupancy(core_result = bcse_core, core_set = "increase")

# Fit and plot neutral model
bcse_neutral <- fit_neutral_model(
  otu_table = bcse_core$otu_table,
  core_set = bcse_core$increase_core,
  abundance_occupancy = bcse_core$abundance_occupancy
)

plot_neutral_model(bcse_neutral)
```

## Documentation

For detailed examples and use cases, see the package vignette:

``` r
vignette("BRCore-vignette", package = "BRCore")
```

## Contributing

Contributions to BRCore are welcome! Please see the
[CONTRIBUTING.md](http://www.germslab.org/BRCore/CONTRIBUTING.md) file
for guidelines on how to contribute.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](http://www.germslab.org/BRCore/CODE_OF_CONDUCT.md). By
participating in this project you agree to abide by its terms.
