# Run multiple rarefaction for microbiome count tables

This function performs multiple rarefaction on a `phyloseq` object by
randomly sub-sampling OTUs/ASVs within samples without replacement. The
process is repeated for a specified number of iterations, and the
results are averaged. Samples with fewer OTUs/ASVs than the specified
`depth_level` are discarded.

## Usage

``` r
multi_rarefy(
  physeq_obj,
  depth_level,
  num_iter = 100,
  .as_array = FALSE,
  set_seed = NULL
)
```

## Arguments

- physeq_obj:

  A `phyloseq` object containing an OTU/ASV table.

- depth_level:

  An integer specifying the sequencing depth (number of OTUs/ASVs) to
  which samples should be rarefied.

- num_iter:

  An integer specifying the number of iterations to perform for
  rarefaction.

- .as_array:

  A logical indicating whether to return the results as a 3D array or as
  a list of data frames. If `TRUE`, returns a 3D array with dimensions
  (samples x taxa x iterations). If `FALSE`, returns a list of data
  frames, one for each iteration, with samples as rows and taxa as
  columns. (default = FALSE)

- set_seed:

  An optional integer to set the random seed for reproducibility
  (default = NULL).

## Value

A data frame with taxa as rows and samples as columns. The values
represent the average sequence counts calculated across all iterations.
Samples with less than `depth_level` sequences are discarded.

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
data("bcse", package = "BRCore")

# Example rarefaction (single iteration, single core to keep examples fast)
otu_table_rare <- multi_rarefy(
  physeq_obj = bcse,
  depth_level = 1000,
  num_iter = 100,
  .as_array = FALSE,
  set_seed = 7642
)
#> Error in multi_rarefy(physeq_obj = bcse, depth_level = 1000, num_iter = 100,     .as_array = FALSE, set_seed = 7642): unused arguments (physeq_obj = bcse, .as_array = FALSE)

rowSums(otu_table_rare[[1]])
#> Error: object 'otu_table_rare' not found
# }
```
