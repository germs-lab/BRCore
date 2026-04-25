# Add a rarefied otu_table to a phyloseq object

This function updates a `phyloseq` object by replacing its OTU/ASV table
with a rarefied version produced by
[`multi_rarefy()`](http://www.germslab.org/BRCore/reference/multi_rarefy.md).
The rarefied table can be a data frame, a list of data frames
(`.as = "list"`), or a 3D array (`.as = "array"`). When providing a list
or array, specify which iteration to use via the `iteration` parameter.

## Usage

``` r
update_otu_table(physeq_obj, rarefied_otus, iteration = NULL)
```

## Arguments

- physeq_obj:

  A `phyloseq` object in which you want to add the rarefied OTU/ASV
  table.

- rarefied_otus:

  A data frame, list of data frames, or 3D array output from
  [`multi_rarefy()`](http://www.germslab.org/BRCore/reference/multi_rarefy.md)
  containing the rarefied OTU/ASV tables.

- iteration:

  Integer specifying which iteration to extract from a list or array.
  Required when `rarefied_otus` is a list or array. Ignored when
  `rarefied_otus` is a data frame.

## Value

A `phyloseq` object.

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
data(GlobalPatterns, package = "phyloseq")

# List output (.as = "list")
otu_list <-
  multi_rarefy(
    physeq_obj = GlobalPatterns,
    depth_level = 200,
    num_iter = 3,
    .as = "list",
    set_seed = 123
  )
#> 
#> ── Multiple Rarefaction ────────────────────────────────────────────────────────
#> 
#> ── Input Validation ──
#> 
#> ✔ Input phyloseq object is valid!
#> ℹ Seed: 123
#> ℹ Input (matrix/df dim): 26 samples x 19216 taxa
#> ℹ Rarefaction depth: 200
#> ℹ Iterations: 3
#> ℹ taxa_are_rows: TRUE
#> ℹ OTU matrix/df rownames head: CL3, CC1, SV1, M31Fcsw, M11Fcsw, M31Plmr
#> ℹ OTU matrix/df colnames head: 549322, 522457, 951, 244423, 586076, 246140
#> ℹ Row sums summary: Min=58688, Max=2357181, Median=1106849
#> 
#> ── Rarefaction Results ──
#> 
#> ── Sample Removal 
#> ✔ No samples removed.
#> 
#> ── Taxa Removal 
#> ✔ No taxa removed.
#> ! Taxa are not removed across iterations to maintain consistent dimensions. 
#> Downstream analyses should handle zero-abundance taxa appropriately.
#> 
#> ── Data Sparsity 
#> ℹ Returning list of data frames for each iteration.
#> • Rarefied matrix (across 3 iterations):
#>   • Min: 497462 zeros (99.57% sparsity) out of 499616 entries
#>   • Max: 497501 zeros (99.58% sparsity) out of 499616 entries
#>   • Avg: 497484.7 zeros (99.57% sparsity) out of 499616 entries
#> 
#> ── Final Data Dimensions 
#> ✔ Output: 3 iterations with 26 unique samples
#> • Samples per iteration:
#>   • Min: 26
#>   • Max: 26
#> • Non-zero taxa per iteration:
#>   • Min: 1476
#>   • Max: 1492
#>   • Avg: 1483.7

# Extract iteration 2
rarefied_gp <- update_otu_table(GlobalPatterns, otu_list, iteration = 2)
#> ℹ Extracting iteration 2 from list of 3 iterations.
#> ✔ Phyloseq object and rarefied otu_table sample names are identical.
#> ✔ All samples kept after rarefaction at depth of: 200
#> ℹ Building phyloseq object with 26 samples and 19216 taxa
#> ✔ Update complete!

# Array output (.as = "array")
otu_array <-
  multi_rarefy(
    physeq_obj = GlobalPatterns,
    depth_level = 200,
    num_iter = 3,
    .as = "array",
    set_seed = 123
  )
#> 
#> ── Multiple Rarefaction ────────────────────────────────────────────────────────
#> 
#> ── Input Validation ──
#> 
#> ✔ Input phyloseq object is valid!
#> ℹ Seed: 123
#> ℹ Input (matrix/df dim): 26 samples x 19216 taxa
#> ℹ Rarefaction depth: 200
#> ℹ Iterations: 3
#> ℹ taxa_are_rows: TRUE
#> ℹ OTU matrix/df rownames head: CL3, CC1, SV1, M31Fcsw, M11Fcsw, M31Plmr
#> ℹ OTU matrix/df colnames head: 549322, 522457, 951, 244423, 586076, 246140
#> ℹ Row sums summary: Min=58688, Max=2357181, Median=1106849
#> 
#> ── Rarefaction Results ──
#> 
#> ── Sample Removal 
#> ✔ No samples removed.
#> 
#> ── Taxa Removal 
#> ✔ No taxa removed.
#> ! Taxa are not removed across iterations to maintain consistent dimensions. 
#> Downstream analyses should handle zero-abundance taxa appropriately.
#> 
#> ── Data Sparsity 
#> • Rarefied matrix (across 3 iterations):
#>   • Min: 497450 zeros (99.57% sparsity) out of 499616 entries
#>   • Max: 497497 zeros (99.58% sparsity) out of 499616 entries
#>   • Avg: 497476.3 zeros (99.57% sparsity) out of 499616 entries
#> 
#> ── Final Data Dimensions 
#> ✔ Output: 3 iterations with 26 samples
#> • Array dimensions:
#>   • Samples: 26
#>   • Taxa: 19216
#>   • Iterations: 3
#> • Non-zero taxa per iteration:
#>   • Min: 1468
#>   • Max: 1511
#>   • Avg: 1483

# Extract iteration 1
rarefied_gp2 <- update_otu_table(GlobalPatterns, otu_array, iteration = 1)
#> ℹ Extracting iteration 1 from array with 3 iterations.
#> ✔ Phyloseq object and rarefied otu_table sample names are identical.
#> ✔ All samples kept after rarefaction at depth of: 200
#> ℹ Building phyloseq object with 26 samples and 19216 taxa
#> ✔ Update complete!
# }
```
