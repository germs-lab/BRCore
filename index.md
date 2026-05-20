# BRCore

## Overview

BRCore provides a unified framework for identification and ecological
interpretation of core microbiomes across time and space, enhancing
robustness and reproducibility in microbiome data analysis.

BRCore provides tools for:

**Rarefaction analysis**: Calculate pre-rarefaction metrics and perform
rarefactions

- [`add_rarefaction_metrics()`](https://www.germslab.org/BRCore/reference/add_rarefaction_metrics.md)
- [`multi_rarefy()`](https://www.germslab.org/BRCore/reference/multi_rarefy.md)
- [`plot_variance_propagation()`](https://www.germslab.org/BRCore/reference/plot_variance_propagation.md)

**Core microbiome identification**: Identify core microbial taxa using
abundance-occupancy distributions

- [`identify_core()`](https://www.germslab.org/BRCore/reference/identify_core.md)

**Neutral model fitting**: Fit and visualize neutral community models

- [`fit_neutral_model()`](https://www.germslab.org/BRCore/reference/fit_neutral_model.md)
- [`plot_neutral_model()`](https://www.germslab.org/BRCore/reference/plot_neutral_model.md)
- [`sncm.fit()`](https://www.germslab.org/BRCore/reference/sncm.fit.md)

**Visualization**: Plot rarefaction diagnostics, abundance-occupancy
curves, and core distributions

- [`plot_rarefaction_metrics()`](https://www.germslab.org/BRCore/reference/plot_rarefaction_metrics.md)
- [`plot_abundance_occupancy()`](https://www.germslab.org/BRCore/reference/plot_abundance_occupancy.md)
- [`plot_core_distribution()`](https://www.germslab.org/BRCore/reference/plot_core_distribution.md)
- [`plot_identified_core()`](https://www.germslab.org/BRCore/reference/plot_identified_core.md)

## Installation

### *Stable* version from CRAN with:

``` r

install.packages("BRCore")

# or if using with `renv`:
renv::init(bioconductor = TRUE)
renv::install("BRCore")
```

***NOTE:*** *If using
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html) in
an `renv` project, `renv` will intercept the call. Run
`renv::init(bioconductor = TRUE)` first to ensure Bioconductor
dependencies can be resolved and installed.*

### *Development* version from GitHub with:

``` r

# install.packages("pak")
pak::pak("germs-lab/BRCore")
```

***Note:*** *`pak` handles dependencies automatically.*

## Quick Start

``` r

library(BRCore)
library(phyloseq)

# Add rarefaction metrics
bcse_metrics <- add_rarefaction_metrics(data = bcse)

# Perform rarefaction
bcse_rarefied_list <- multi_rarefy(
  physeq_obj = bcse,
  depth_level = 1000,
  num_iter = 3,
  set_seed = 7642
)

# Update phyloseq object with rarefied data
bcse_rare_single <- update_otu_table(
  physeq_obj = bcse,
  rarefied_otus = bcse_rarefied_list,
  iteration = 2
) # Your preffered iteration can be used here

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
  rarefied_list = bcse_rarefied_list,
  priority_var = "Crop",
  increase_value = 0.02,
  depth_level = 1000,
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
[CONTRIBUTING.md](https://github.com/germs-lab/BRCore/blob/main/.github/CONTRIBUTING.md)
file for guidelines on how to contribute.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/germs-lab/BRCore/blob/main/.github/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
