# Identify Core Microbial Taxa

This function identifies core microbial taxa based on
abundance-occupancy distributions and their contributions to Bray-Curtis
similarity between biological samples. Core taxa are selected using
either a "last % increase" or "elbow" method implementing the method
developed by Shade and Stopnisek (2019) Curr Opin Microbiol, see below
for details.

## Usage

``` r
identify_core(
  physeq_obj,
  priority_var,
  increase_value = 0.02,
  abundance_weight = 0,
  max_otus = NULL,
  depth_level = 1000,
  num_iter = 100,
  seed = NULL
)
```

## Arguments

- physeq_obj:

  A `phyloseq` object with at least `otu_table` and `sample_data`.

- priority_var:

  The column name in the `sample_data` (e.g. sampling_date", "site")
  that is used for prioritizing the core microbiome.

- increase_value:

  Increase value (numeric, scalar) used in the calculation (default
  0.02) for "increase". The "elbow" is always calculated and returned as
  `elbow_core` (see below for details).

- abundance_weight:

  Numeric in `[0,1]`; how much to weight mean relative abundance in the
  ranking score. `0` (default) uses occupancy/composite only. `1` ranks
  purely by abundance. Values in between blend the two (e.g.,
  abundance_weight = 0.3 gives 70% occupancy/composite + 30% abundance).

- max_otus:

  Optional integer to limit analysis to the top N ranked OTUs. If NULL
  (default), all OTUs are analyzed. Useful for large datasets (\>5000
  OTUs)

- depth_level:

  Integer. The sequencing depth used for normalization in Bray-Curtis
  calculations. If data is rarefied, this is automatically set to the
  rarefaction depth. For unrarefied data, samples with depth below this
  threshold are excluded from pairwise comparisons.

- num_iter:

  Integer. Number of subsampling iterations used when calculating
  average dissimilarity for unrarefied data (default 100). Ignored if
  data is already rarefied.

- seed:

  Optional integer to set the RNG seed for reproducibility.

## Value

A list with:

- `bray_curtis_ranked` tibble with `rank`, mean Bray-Curtis similarity
  across sample pairs (`MeanBC`) at each cumulative `rank`, normalized
  proportion (`proportionBC`), the multiplicative `IncreaseBC`, and the
  elbow metric (`elbow_slope_diffs`). (`proportionBC`), the
  multiplicative `IncreaseBC`, and the elbow metric
  (`elbow_slope_diffs`).

- `otu_ranked` tibble with ranked OTU/ASVs .

- `abundance_occupancy` tibble with OTU/ASVs names, occupancy
  (`otu_occ`), and mean relative abundance (`otu_rel`).

- `priority_var` character, the variable used for prioritizing the core.

- `elbow` core set identified by elbow method (integer).

- `bc_increase` core set identified by last % BC-increase (integer).

- `increase_value` increase value (numeric, scalar) used in the
  calculation (e.g. 0.02).

- `elbow_core` core OTU/ASVs using elbow method (character vector).

- `increase_core` core OTU/ASVs using last % BC-increase method
  (character vector).

- `otu_table` otu_table counts (otu x samples) used (data.frame).

- `sample_metadata` samples metadata (data.frame).

- `taxonomy_table` taxonomy if present (data.frame); otherwise NULL.

## Details

The core set is defined using two separate methods:

The function rank OTU/ASVs by occupancy (optionally with abundance
weighting:
`rank_score = (1 - weight) * Index + weight * scaled_abundance`, where
scaled_abundance is mean relative abundance rescaled to `[0,1]`). For
each `k = 1..K`, recompute `S_k` as the mean Bray-Curtis similarity
across all sample pairs using only the first `k` ranked OTUs; when
`k = K`, this yields `S_K`, the value computed with all OTUs.
Normalizing by `S_K` gives `C_k = S_k / S_K`.

The **elbow** is the point of diminishing returns: for each `k`, compare
the average *left* slope `(S_k - S_1) / (k - 1)` to the average *right*
slope `(S_K - S_k) / (K - k)`, and choose the `k` that maximizes
`(left - right)`.

The **last percent Bray-Curtis increase** method uses the same
accumulation curve, examine the multiplicative step when adding the
`k-th` OTU: `Increase_k = S_k / S_{k-1}` (equivalently,
`Increase_k = C_k / C_{k-1}`). Choose the largest `k` such that
`Increase_k >= 1 + p`, where `p` is your chosen percent threshold
(increase_value; recommended `p >= 0.02` or `2%`). This selects the
final rank for which adding one more taxon still increases the explained
Bray-Curtis similarity by at least `p`.

## Dependencies

Requires phyloseq, dplyr, tidyr, tibble, rlang, and vegan.

## References

Shade A, Stopnisek N (2019) Abundance-occupancy distributions to
prioritize plant core microbiome membership. Current Opinion in
Microbiology, 49:50-58 doi:https://doi.org/10.1016/j.mib.2019.09.008

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
# Example using your switchgrass phyloseq object and grouping variable
# 'sampling_date'
data("switchgrass", package = "BRCore")

res <- identify_core(
  physeq_obj     = switchgrass,
  priority_var   = "sampling_date",
  increase_value = 0.02,
  seed           = 091825
)
#> Seed used: 91825
#> ✔ Input phyloseq object is valid!
#> ℹ otu_table() is rarefied at a depth of: 1000
#> ℹ No taxonomy found (or empty). Continuing without taxonomy.
#> ✔ Core prioritizing variable: sampling_date
#> ℹ Ranking OTUs based on BC dissimilarity, starting at 2026-04-21 12:00:23.380814
#> ✔ Elbow method identified 3 core OTUs
#> ✔ % increase method identified 34 core OTUs
#> ✔ Analysis complete!

# Inspect results
str(res)
#> List of 13
#>  $ bray_curtis_ranked : tibble [706 × 5] (S3: tbl_df/tbl/data.frame)
#>   ..$ rank             : Factor w/ 706 levels "1","10","100",..: 1 112 223 334 445 556 667 685 696 2 ...
#>   ..$ MeanBC           : num [1:706] 0.00638 0.04689 0.09625 0.10422 0.10888 ...
#>   ..$ proportionBC     : num [1:706] 0.0124 0.0913 0.1875 0.203 0.2121 ...
#>   ..$ IncreaseBC       : num [1:706] 0 7.35 2.05 1.08 1.04 ...
#>   ..$ elbow_slope_diffs: num [1:706] -0.000719 0.019593 0.029361 0.023876 0.019922 ...
#>  $ otu_ranked         :'data.frame': 706 obs. of  8 variables:
#>   ..$ otu    : chr [1:706] "OTU47" "OTU2" "OTU6" "OTU21" ...
#>   ..$ otu_occ: num [1:706] 1 1 1 1 1 ...
#>   ..$ otu_rel: num [1:706] 0.0161 0.1677 0.165 0.0131 0.0114 ...
#>   ..$ sumF   : num [1:706] 6 6 6 6 6 ...
#>   ..$ sumG   : num [1:706] 6 6 6 6 6 5 5 5 5 5 ...
#>   ..$ nS     : int [1:706] 12 12 12 12 12 12 12 12 12 12 ...
#>   ..$ Index  : num [1:706] 1 1 1 1 1 ...
#>   ..$ rank   : num [1:706] 1 1 1 1 1 ...
#>  $ abundance_occupancy:'data.frame': 706 obs. of  3 variables:
#>   ..$ otu    : chr [1:706] "OTU10713" "OTU7" "OTU1253" "OTU47" ...
#>   ..$ otu_occ: num [1:706] 0.2093 0.9767 0.0698 1 0.9535 ...
#>   ..$ otu_rel: num [1:706] 0.000233 0.046767 0.000419 0.016093 0.011953 ...
#>  $ priority_var       : chr "sampling_date"
#>  $ abundance_weight   : num 0
#>  $ elbow              : int 3
#>  $ bc_increase        : int 34
#>  $ increase_value     : num 0.02
#>  $ elbow_core         : chr [1:3] "OTU47" "OTU2" "OTU6"
#>  $ increase_core      : chr [1:34] "OTU47" "OTU2" "OTU6" "OTU21" ...
#>  $ otu_table          : int [1:706, 1:43] 0 129 0 10 23 1 0 1 33 276 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:706] "OTU10713" "OTU7" "OTU1253" "OTU47" ...
#>   .. ..$ : chr [1:43] "G5R1_MAIN_01AUG2016_LD2" "G5R1_MAIN_03OCT2016_LD1" "G5R1_MAIN_12JUL2016_LD1" "G5R1_MAIN_12SEP2016_LD2" ...
#>  $ metadata           :'data.frame': 43 obs. of  3 variables:
#>   ..$ plant        : chr [1:43] "switchgrass" "switchgrass" "switchgrass" "switchgrass" ...
#>   ..$ sampling_date: POSIXct[1:43], format: "2016-07-31 23:00:00" "2016-10-02 23:00:00" ...
#>   ..$ sample_id    : chr [1:43] "G5R1_MAIN_01AUG2016_LD2" "G5R1_MAIN_03OCT2016_LD1" "G5R1_MAIN_12JUL2016_LD1" "G5R1_MAIN_12SEP2016_LD2" ...
#>  $ taxonomy           : NULL
#>  - attr(*, "class")= chr [1:2] "identify_core_result" "list"
# }
```
