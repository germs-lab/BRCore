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

## Examples

``` r
# \donttest{
library(phyloseq)
library(BRCore)
# From an object class "phyloseq" with added alpha metrics

data("bcse", package = "BRCore")
bcse_metrics <- add_rarefaction_metrics(data = bcse)
sample_data(bcse_metrics)
#>         Niche               Crop Plot read_num singlton_num goods_cov outlier
#> bcse50   Leaf               Corn   R2   107643          160  99.85136      NA
#> bcse69   Leaf            Sorghum   R4    54331          190  99.65029      NA
#> bcse73   Leaf        Switchgrass   R5    67815           82  99.87908      NA
#> bcse191  Leaf         Miscanthus   R5    83152           96  99.88455      NA
#> bcse82   Leaf     Native Grasses   R2    32228          156  99.51595      NA
#> bcse102  Leaf             Poplar   R3     9621           52  99.45952      NA
#> bcse111  Leaf   Early Succession   R4     6699           43  99.35811      NA
#> bcse86   Leaf         Miscanthus   R3   100940           85  99.91579      NA
#> bcse97   Leaf   Early Succession   R3     4020           30  99.25373      NA
#> bcse88   Leaf         Miscanthus   R1    70532           90  99.87240      NA
#> bcse104  Leaf             Poplar   R1     7819           28  99.64190      NA
#> bcse108  Leaf            Prairie   R2      627           29  95.37480      NA
#> bcse81   Leaf     Native Grasses   R3     5532           66  98.80694      NA
#> bcse77   Leaf        Switchgrass   R2     6425           45  99.29961      NA
#> bcse78   Leaf        Switchgrass   R1    22305           93  99.58305      NA
#> bcse96   Leaf     Native Grasses   R4    19814           83  99.58110      NA
#> bcse66   Leaf    New Switchgrass   R2    18536          305  98.35455      NA
#> bcse57   Leaf Continuous Sorghum   R5    60437          169  99.72037      NA
#> bcse75   Leaf        Switchgrass   R4     3146           59  98.12460      NA
#> bcse101  Leaf             Poplar   R4     4529           27  99.40384      NA
#> bcse192  Leaf               Corn   R4    89228          161  99.81956      NA
#> bcse98   Leaf   Early Succession   R2     2649           40  98.49000      NA
#> bcse105  Leaf            Prairie   R5      151           19  87.41722      NA
#> bcse106  Leaf            Prairie   R4     1354           32  97.63663      NA
#> bcse76   Leaf        Switchgrass   R3     4936           33  99.33144      NA
#> bcse103  Leaf             Poplar   R2     9552           39  99.59171      NA
#> bcse51   Leaf               Corn   R1    94449          216  99.77131      NA
#> bcse63   Leaf               Corn   R5    91353          141  99.84565      NA
#> bcse68   Leaf            Sorghum   R5    78348          120  99.84684      NA
#> bcse109  Leaf            Prairie   R1     3229           80  97.52245      NA
#> bcse65   Leaf    New Switchgrass   R3     1193          165  86.16932      NA
#> bcse58   Leaf Continuous Sorghum   R4    53594          102  99.80968      NA
#> bcse107  Leaf            Prairie   R3     1725          313  81.85507      NA
#> bcse62   Leaf Continuous Sorghum   R1    79380           78  99.90174      NA
#> bcse59   Leaf Continuous Sorghum   R3    90635          261  99.71203      NA
#> bcse99   Leaf   Early Succession   R1    60638           59  99.90270      NA
#> bcse49   Leaf               Corn   R3    88058          109  99.87622      NA
#> bcse85   Leaf         Miscanthus   R4    75760           64  99.91552      NA
#> bcse79   Leaf    New Switchgrass   R5     4275          224  94.76023      NA
#> bcse72   Leaf            Sorghum   R1    76278          138  99.81908      NA
#> bcse80   Leaf    New Switchgrass   R4     1963          137  93.02089      NA
#> bcse71   Leaf            Sorghum   R2   102787           77  99.92509      NA
#> bcse110  Leaf   Early Succession   R5      944           27  97.13983      NA
#> bcse95   Leaf     Native Grasses   R5    21592           61  99.71749      NA
#> bcse67   Leaf    New Switchgrass   R1    24099          211  99.12444      NA
#> bcse61   Leaf Continuous Sorghum   R2    86626          114  99.86840      NA
#> bcse100  Leaf             Poplar   R5    16463           44  99.73273      NA
#> bcse87   Leaf         Miscanthus   R2    74457           56  99.92479      NA
#> bcse83   Leaf     Native Grasses   R1    25209          148  99.41291      NA
#> bcse70   Leaf            Sorghum   R3    66743          145  99.78275      NA

# From a class "data.frame" count table object

bcse_otutable <- as.data.frame(
  as(otu_table(bcse), "matrix")
)
test_otutable_metrics <- add_rarefaction_metrics(
  data = bcse_otutable
)
test_otutable_metrics[
  utils::tail(seq_len(nrow(test_otutable_metrics)), 10),
  utils::tail(seq_len(ncol(test_otutable_metrics)), 20)
]
#>           bcse59 bcse99 bcse49 bcse85 bcse79 bcse72 bcse80 bcse71 bcse110
#> OTU_15053      1      0      0      0      0      0      0      0       0
#> OTU_14114      0      0      0      0      0      0      0      0       0
#> OTU_9931       0      0      0      0      0      0      0      0       0
#> OTU_13974      2      0      0      0      0      0      0      0       0
#> OTU_10518      0      0      0      0      0      0      0      6       0
#> OTU_12854      0      0      0      0      0      0      0      1       0
#> OTU_15056      0      0      0      0      0      2      0      0       0
#> OTU_14409      0      0      0      0      0      0      0      0       0
#> OTU_11281      0      0      0      0      0      0      0      0       0
#> OTU_15187      0      0      0      2      0      0      0      0       0
#>           bcse95 bcse67 bcse61 bcse100 bcse87 bcse83 bcse70 read_num
#> OTU_15053      0      0      0       0      0      0      0        2
#> OTU_14114      0      0      0       0      0      0      0        1
#> OTU_9931       0      0      0       0      0      0      0        1
#> OTU_13974      0      0      0       0      0      0      0        2
#> OTU_10518      0      0      0       0      0      0      0        6
#> OTU_12854      0      0      0       0      0      0      0        2
#> OTU_15056      0      0      0       0      0      0      0        2
#> OTU_14409      0      0      0       0      0      0      0        2
#> OTU_11281      0      0      0       0      0      0      0        2
#> OTU_15187      0      0      0       0      0      0      0        2
#>           singlton_num goods_cov outlier
#> OTU_15053            2         0      NA
#> OTU_14114            1         0      NA
#> OTU_9931             1         0      NA
#> OTU_13974            0       100      NA
#> OTU_10518            0       100      NA
#> OTU_12854            2         0      NA
#> OTU_15056            0       100      NA
#> OTU_14409            0       100      NA
#> OTU_11281            0       100      NA
#> OTU_15187            0       100      NA
# }
```
