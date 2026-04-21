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
#> 
#> ── Multiple Rarefaction ────────────────────────────────────────────────────────
#> 
#> ── Input Validation ──
#> 
#> ✔ Input phyloseq object is valid!
#> ℹ Seed: 7642
#> ℹ Input (matrix/df dim): 47 samples x 2861 taxa
#> ℹ Rarefaction depth: 1000
#> ℹ Iterations: 100
#> ℹ taxa_are_rows: TRUE
#> ℹ OTU matrix/df rownames head: bcse50, bcse69, bcse73, bcse191, bcse82, bcse102
#> ℹ OTU matrix/df colnames head: OTU_427, OTU_11, OTU_253, OTU_148, OTU_3, OTU_78
#> ℹ Row sums summary: Min=1193, Max=107643, Median=25209
#> 
#> ── Rarefaction Results ──
#> 
#> ── Sample Removal 
#> ! 3 samples removed (depth < 1000)
#> ! Samples removed: "bcse108, bcse105, bcse110"
#> 
#> ── Taxa Removal 
#> ℹ Original taxa input: 2861
#> ! Max: 1832 taxa removed (zero abundance) in viable samples (depth_level >= 1000).
#> ! When using `.as_array = FALSE`, taxa removed may differ across iterations.
#> 
#> ── Data Sparsity 
#> ℹ Returning list of data frames for each iteration.
#> • Rarefied matrix (across 100 iterations):
#>   • Min: 44542 zeros (92.1% sparsity) out of 48363 entries
#>   • Max: 48865 zeros (92.58% sparsity) out of 52781 entries
#>   • Avg: 46569.7 zeros (92.36% sparsity) out of 50420.7 entries
#> 
#> ── Final Data Dimensions 
#> ✔ Output: 100 iterations with 50 unique samples
#> • Samples per iteration:
#>   • Min: 47
#>   • Max: 47
#> • Non-zero taxa per iteration:
#>   • Min: 1029
#>   • Max: 1123
#>   • Avg: 1072.8

rowSums(otu_table_rare[[1]])
#>  bcse50  bcse69  bcse73 bcse191  bcse82 bcse102 bcse111  bcse86  bcse97  bcse88 
#>    1000    1000    1000    1000    1000    1000    1000    1000    1000    1000 
#> bcse104  bcse81  bcse77  bcse78  bcse96  bcse66  bcse57  bcse75 bcse101 bcse192 
#>    1000    1000    1000    1000    1000    1000    1000    1000    1000    1000 
#>  bcse98 bcse106  bcse76 bcse103  bcse51  bcse63  bcse68 bcse109  bcse65  bcse58 
#>    1000    1000    1000    1000    1000    1000    1000    1000    1000    1000 
#> bcse107  bcse62  bcse59  bcse99  bcse49  bcse85  bcse79  bcse72  bcse80  bcse71 
#>    1000    1000    1000    1000    1000    1000    1000    1000    1000    1000 
#>  bcse95  bcse67  bcse61 bcse100  bcse87  bcse83  bcse70 
#>    1000    1000    1000    1000    1000    1000    1000 
# }
```
