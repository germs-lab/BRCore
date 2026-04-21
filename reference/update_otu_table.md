# Add a rarefied otu_table to a phyloseq object

This function updates a `phyloseq` object by replacing its OTU/ASV table
with a rarefied version produced by
[`multi_rarefy()`](http://www.germslab.org/BRCore/reference/multi_rarefy.md).
The rarefied table can be a data frame, a list of data frames
(`.as_array = FALSE`), or a 3D array (`.as_array = TRUE`). When
providing a list or array, specify which iteration to use via the
`iteration` parameter.

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

# List output (.as_array = FALSE)
otu_list <-
  multi_rarefy(
    physeq_obj = GlobalPatterns,
    depth_level = 200,
    num_iter = 3,
    .as_array = FALSE,
    set_seed = 123
  )
#> Error in multi_rarefy(physeq_obj = GlobalPatterns, depth_level = 200,     num_iter = 3, .as_array = FALSE, set_seed = 123): unused arguments (physeq_obj = GlobalPatterns, .as_array = FALSE)

# Extract iteration 2
rarefied_gp <- update_otu_table(GlobalPatterns, otu_list, iteration = 2)
#> Error in update_otu_table(GlobalPatterns, otu_list, iteration = 2): unused argument (iteration = 2)

# Array output (.as_array = TRUE)
otu_array <-
  multi_rarefy(
    physeq_obj = GlobalPatterns,
    depth_level = 200,
    num_iter = 3,
    .as_array = TRUE,
    set_seed = 123
  )
#> Error in multi_rarefy(physeq_obj = GlobalPatterns, depth_level = 200,     num_iter = 3, .as_array = TRUE, set_seed = 123): unused arguments (physeq_obj = GlobalPatterns, .as_array = TRUE)

# Extract iteration 1
rarefied_gp2 <- update_otu_table(GlobalPatterns, otu_array, iteration = 1)
#> Error in update_otu_table(GlobalPatterns, otu_array, iteration = 1): unused argument (iteration = 1)
# }
```
