# Calculate and append pre-rarefaction statistics to microbiome data

This function adds read count, singleton count, Good's coverage, and
marks outlier samples to a `phyloseq` object or `data.frame` based on
the OTU/ASV abundance table.

## Usage

``` r
add_rarefaction_metrics(data)
```

## Arguments

- data:

  A `phyloseq` object or a `data.frame` with samples as rows and taxa as
  columns.

## Value

The same object (`phyloseq` or `data.frame`) with new columns:

- `read_num`

- `singleton_num`

- `goods_cov`

- `outlier`

## Details

About Good's coverage. Initially developed by Alan Turing and I.J. Good
during their cryptographic analyses in World War II, it was later
adopted by ecologists, particularly in microbial diversity studies, to
assess the completeness of a sample's representation of the overall
community. It's calculated as `1 - (F1/N)`, where `F1` is the number of
OTUs (Operational Taxonomic Units) represented by only one individual
(singletons) and `N` is the total number of individuals in the sample.
For example, a Good's coverage of 0.95, means that 5% of the reads in
that sample are from OTUs that appear only once.

## See also

[`plot_rarefaction_metrics()`](http://www.germslab.org/BRCore/reference/plot_rarefaction_metrics.md)
for visualizing these metrics, and
[`multi_rarefy()`](http://www.germslab.org/BRCore/reference/multi_rarefy.md)
for performing rarefaction on a `phyloseq` object.

## Examples

``` r
library(phyloseq)
library(BRCore)

data("bcse", package = "BRCore")

# Adding metrics to a "phyloseq" object
bcse_metrics <- add_rarefaction_metrics(data = bcse)
sample_data(bcse_metrics)|>
head(10)
#>         Niche             Crop Plot read_num singleton_num goods_cov outlier
#> bcse50   Leaf             Corn   R2   107643           160  99.85136      NA
#> bcse69   Leaf          Sorghum   R4    54331           190  99.65029      NA
#> bcse73   Leaf      Switchgrass   R5    67815            82  99.87908      NA
#> bcse191  Leaf       Miscanthus   R5    83152            96  99.88455      NA
#> bcse82   Leaf   Native Grasses   R2    32228           156  99.51595      NA
#> bcse102  Leaf           Poplar   R3     9621            52  99.45952      NA
#> bcse111  Leaf Early Succession   R4     6699            43  99.35811      NA
#> bcse86   Leaf       Miscanthus   R3   100940            85  99.91579      NA
#> bcse97   Leaf Early Succession   R3     4020            30  99.25373      NA
#> bcse88   Leaf       Miscanthus   R1    70532            90  99.87240      NA


# Adding metrics to a "data.frame" count table object
bcse_otutable <- as.data.frame(
  as.matrix(otu_table(bcse))
)

bcse_otutable_metrics <- add_rarefaction_metrics(
  data = bcse_otutable
)
bcse_otutable_metrics[
  head(seq_len(nrow(bcse_otutable_metrics)), 10),
  tail(seq_len(ncol(bcse_otutable_metrics)), 20)
]
#>         bcse59 bcse99 bcse49 bcse85 bcse79 bcse72 bcse80 bcse71 bcse110 bcse95
#> OTU_427      2      0      0      0      0      0      1      0       0      0
#> OTU_11       5      0      0      0     12      0     37      0       1      2
#> OTU_253      0      0      0      0      0      0      0      0       0      0
#> OTU_148    323      0      6      1     15    711      7   1107       0      0
#> OTU_3    30332  11876  34910  23513    111   6016    201   1130     163   3992
#> OTU_78       0      1      0      0      1      2      0      0       0      0
#> OTU_58       9      0      0      0      6      0      1      0       0      1
#> OTU_152      0      0      0      0      0      0      0      5       0      0
#> OTU_16      37    195     30     12     46     24      3      3       3     31
#> OTU_35       0      0      0      0      0      0      0      0       0      0
#>         bcse67 bcse61 bcse100 bcse87 bcse83 bcse70 read_num singleton_num
#> OTU_427      3      0       0      0      0      0        7             2
#> OTU_11       3      1       0      0     29      0      127            12
#> OTU_253      2      0       0      0      0      0        5             0
#> OTU_148     13   5093     184     34     70     13    16078             9
#> OTU_3     8004  26526    3140  19678   5773   2838   510842             0
#> OTU_78       1      0       0      0      0      0       30            10
#> OTU_58      16      0       0      0      0      3       79             5
#> OTU_152      0      0       0      0      0      0        8             1
#> OTU_16      58     64      19      2     95     19     1463             5
#> OTU_35       0      0       0      0      0      1        2             2
#>         goods_cov outlier
#> OTU_427  71.42857      NA
#> OTU_11   90.55118      NA
#> OTU_253 100.00000      NA
#> OTU_148  99.94402   16078
#> OTU_3   100.00000  510842
#> OTU_78   66.66667      NA
#> OTU_58   93.67089      NA
#> OTU_152  87.50000      NA
#> OTU_16   99.65824    1463
#> OTU_35    0.00000      NA

```
