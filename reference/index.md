# Package index

## Core Identification Functions

- [`identify_core()`](https://www.germslab.org/BRCore/reference/identify_core.md)
  : Identify Core Microbiome Using Bray-Curtis Similarity biological
  samples. Core taxa are selected using either a "last % increase" or
  "elbow" method implementing the method developed by Shade and
  Stopnisek (2019) Curr Opin Microbiol, see below for details.
- [`plot_identified_core()`](https://www.germslab.org/BRCore/reference/plot_identified_core.md)
  : Plot Bray-Curtis increase over ranked OTU/ASVs
- [`plot_core_distribution()`](https://www.germslab.org/BRCore/reference/plot_core_distribution.md)
  : Plot Core Taxa Occupancy Across Metadata Groups

## Data Exploration Functions

- [`plot_rarefaction_metrics()`](https://www.germslab.org/BRCore/reference/plot_rarefaction_metrics.md)
  : Plot pre-rarefaction diagnostics
- [`plot_abundance_occupancy()`](https://www.germslab.org/BRCore/reference/plot_abundance_occupancy.md)
  : Plot Abundance-Occupancy Curve and Display the Core Taxa
- [`plot_core_distribution()`](https://www.germslab.org/BRCore/reference/plot_core_distribution.md)
  : Plot Core Taxa Occupancy Across Metadata Groups
- [`plot_identified_core()`](https://www.germslab.org/BRCore/reference/plot_identified_core.md)
  : Plot Bray-Curtis increase over ranked OTU/ASVs
- [`plot_variance_propagation()`](https://www.germslab.org/BRCore/reference/plot_variance_propagation.md)
  : Variance propagation diagnostic for rarefaction
- [`plot_neutral_model()`](https://www.germslab.org/BRCore/reference/plot_neutral_model.md)
  : Plot a fitted Neutral Model to Microbial Community Data

## Neutral Model Functions

- [`fit_neutral_model()`](https://www.germslab.org/BRCore/reference/fit_neutral_model.md)
  : Fit a Neutral Model to Microbial Community Data
- [`plot_neutral_model()`](https://www.germslab.org/BRCore/reference/plot_neutral_model.md)
  : Plot a fitted Neutral Model to Microbial Community Data
- [`sncm.fit()`](https://www.germslab.org/BRCore/reference/sncm.fit.md)
  : Fit Sloan Neutral Community Model (SNCM)

## Rarefaction Functions

- [`add_rarefaction_metrics()`](https://www.germslab.org/BRCore/reference/add_rarefaction_metrics.md)
  : Calculate and append pre-rarefaction statistics to microbiome data
- [`multi_rarefy()`](https://www.germslab.org/BRCore/reference/multi_rarefy.md)
  : Run rarefaction for microbiome count tables
- [`plot_rarefaction_metrics()`](https://www.germslab.org/BRCore/reference/plot_rarefaction_metrics.md)
  : Plot pre-rarefaction diagnostics

## OTU Table Utilities

- [`update_otu_table()`](https://www.germslab.org/BRCore/reference/update_otu_table.md)
  : Add a rarefied otu_table to a phyloseq object

## Datasets

- [`bcse`](https://www.germslab.org/BRCore/reference/bcse.md) : 16S
  amplicon dataset from the GLBRC Biofuel Cropping System Experiment
  (BCSE) at Michigan State University, Kellogg Biological Station
- [`bean`](https://www.germslab.org/BRCore/reference/bean.md) : 16S
  amplicon dataset from common bean (Phaseolus vulgaris)
- [`mimulus`](https://www.germslab.org/BRCore/reference/mimulus.md) :
  16S amplicon dataset from yellow monkeyflower (Mimulus guttatus)
- [`switchgrass`](https://www.germslab.org/BRCore/reference/switchgrass.md)
  : 16S amplicon dataset from Switchgrass (Panicum virgatum)
- [`switchgrass_core`](https://www.germslab.org/BRCore/reference/switchgrass_core.md)
  : Identified core microbiome members for the switchgrass dataset
