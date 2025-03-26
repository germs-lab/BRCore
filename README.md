
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRCore

<!-- badges: start -->

[![R-CMD-check](https://github.com/germs-lab/BRCore/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/germs-lab/BRCore/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goal of BRCore is to aid the analysis of the [Inter BRC Microbiome
project](https://github.com/germs-lab/interbrc-core-analysis) by
providing a set of tools to process and analyze microbial data from
Bioenergy Research Centers. ***This is still an experimental package***.

## Installation

Install the development version of `BRCore` from GitHub with:

``` r
# install.packages("pak")
pak::pak("germs-lab/BRCore")

#> ℹ Loading metadata database✔ Loading metadata database ... done
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + BRCore   0.0.0.9000 [bld][cmp] (GitHub: 1ad972a)
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got BRCore 0.0.0.9000 (source) (587.62 kB)
#> ℹ Packaging BRCore 0.0.0.9000
#> ✔ Packaged BRCore 0.0.0.9000 (475ms)
#> ℹ Building BRCore 0.0.0.9000
#> ✔ Built BRCore 0.0.0.9000 (9.5s)
#> ✔ Installed BRCore 0.0.0.9000 (github::germs-lab/BRCore@1ad972a) (25ms)
#> ✔ 1 pkg + 114 deps: kept 101, added 1, dld 1 (NA B) [15.8s]

```

## Step 1: Load Required Libraries and Data

``` r
library(phyloseq)
library(BRCore)

# Load the esophagus dataset
data(esophagus, package = "phyloseq")
```

## Step 2: Create Phyloseq Object

``` r
# Get the taxa names from the esophagus dataset
taxa_names <- taxa_names(esophagus)

# Define realistic taxonomy levels
  kingdoms <- c("Bacteria", "Archaea")
  phyla <- c(
    "Firmicutes",
    "Bacteroidetes",
    "Proteobacteria",
    "Actinobacteria",
    "Euryarchaeota"
  )
  classes <- c(
    "Clostridia",
    "Bacteroidia",
    "Gammaproteobacteria",
    "Actinobacteria",
    "Methanobacteria"
  )
  orders <- c(
    "Clostridiales",
    "Bacteroidales",
    "Enterobacterales",
    "Bifidobacteriales",
    "Methanobacteriales"
  )
  families <- c(
    "Lachnospiraceae",
    "Bacteroidaceae",
    "Enterobacteriaceae",
    "Bifidobacteriaceae",
    "Methanobacteriaceae"
  )
  genera <- c("Blautia",
              "Bacteroides",
              "Escherichia",
              "Bifidobacterium",
              "Methanobrevibacter")
  species <- c(
    "Blautia producta",
    "Bacteroides fragilis",
    "Escherichia coli",
    "Bifidobacterium longum",
    "Methanobrevibacter smithii"
  )
  
  # Create a mock taxonomy table
  set.seed(8998); mock_taxonomy_table <- matrix(
    c(
      sample(kingdoms, length(taxa_names), replace = TRUE),
      sample(phyla, length(taxa_names), replace = TRUE),
      sample(classes, length(taxa_names), replace = TRUE),
      sample(orders, length(taxa_names), replace = TRUE),
      sample(families, length(taxa_names), replace = TRUE),
      sample(genera, length(taxa_names), replace = TRUE),
      sample(species, length(taxa_names), replace = TRUE)
    ),
    nrow = length(taxa_names),
    ncol = 7,
    dimnames = list(
      taxa_names,
      c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      )
    )
  )

# Convert to a taxonomy table object
tax_table <- tax_table(mock_taxonomy_table)

# Add sample metadata to the esophagus dataset
sample_data <- data.frame(
  Sample = sample_names(esophagus),
  Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE), # Random groups
  row.names = sample_names(esophagus)
)

# Add the taxonomy table and sample metadata to the esophagus dataset
esophagus_with_tax <- merge_phyloseq(esophagus, tax_table, sample_data(sample_data))
```

## Step 2: Run extract_core()

``` r
# Extract core microbial taxa
core_result <- extract_core(
  physeq = esophagus_with_tax,
  Var = "Group",  # Variable in sample metadata to group by
  method = "increase",  # Method for core taxa selection
  increase_value = 2  # Threshold for core taxa inclusion (2% contribution to Bray-Curtis dissimilarity)
)
#> ✔ Input phyloseq object is valid!
#> ℹ Calculating BC dissimilarity based on the 1st ranked OTU
#> ✔ BC dissimilarity based on the 1st ranked OTU complete
#> ℹ Calculating BC dissimilarity based on ranked OTUs, starting at 2025-03-26 11:12:55.287585
#> ✔ BC ranks done!
#> ✔ Performing method 'increase'

# View the results
print(core_result)
#> $core_otus
#>  [1] "65_7_12" "65_2_17" "65_7_4"  "59_5_2"  "65_2_5"  "59_7_6"  "59_9_31"
#>  [8] "59_5_19" "59_2_6"  "65_6_2"  "65_3_18" "59_9_17" "59_8_12" "59_4_25"
#> [15] "59_3_21"
#> 
#> $bray_curtis_ranked
#> # A tibble: 44 × 4
#>    rank  MeanBC proportionBC IncreaseBC
#>    <fct>  <dbl>        <dbl>      <dbl>
#>  1 2     0.0985        0.186       1.15
#>  2 3     0.108         0.204       1.1 
#>  3 4     0.117         0.220       1.08
#>  4 5     0.128         0.241       1.10
#>  5 6     0.133         0.251       1.04
#>  6 7     0.169         0.319       1.27
#>  7 8     0.181         0.341       1.07
#>  8 9     0.213         0.402       1.18
#>  9 10    0.351         0.663       1.65
#> 10 11    0.363         0.684       1.03
#> # ℹ 34 more rows
#> 
#> $otu_rankings
#>        otu      rank
#> 1  59_8_22 1.0000000
#> 2  65_7_12 1.0000000
#> 3  65_2_17 1.0000000
#> 4   65_7_4 1.0000000
#> 5   59_5_2 1.0000000
#> 6   65_2_5 1.0000000
#> 7   59_7_6 1.0000000
#> 8  59_9_31 1.0000000
#> 9  59_5_19 1.0000000
#> 10  59_2_6 1.0000000
#> 11  65_6_2 0.3333333
#> 12  9_7_25 0.3333333
#> 13 65_9_13 0.3333333
#> 14 65_3_18 0.3333333
#> 15 59_9_17 0.3333333
#> 16   9_4_6 0.3333333
#> 17 59_5_13 0.1666667
#> 18 59_8_12 0.1666667
#> 19 65_3_22 0.1666667
#> 20  65_5_1 0.1666667
#> 21 65_1_10 0.1666667
#> 22  59_6_1 0.1666667
#> 23 65_9_26 0.1666667
#> 24 65_5_18 0.1666667
#> 25 65_8_12 0.1666667
#> 26  65_9_1 0.1666667
#> 27 59_9_26 0.1666667
#> 28  9_6_28 0.1666667
#> 29 65_4_26 0.1666667
#> 30 65_1_17 0.1666667
#> 31   9_4_3 0.1666667
#> 32  65_8_7 0.1666667
#> 33  59_4_5 0.1666667
#> 34 65_4_10 0.1666667
#> 35   9_4_5 0.1666667
#> 36  65_1_8 0.1666667
#> 37  9_4_13 0.1666667
#> 38   9_6_3 0.1666667
#> 39 59_9_18 0.1666667
#> 40 65_8_25 0.1666667
#> 41 65_8_29 0.1666667
#> 42 65_4_20 0.1666667
#> 43 59_4_25 0.1666667
#> 44 59_3_21 0.1666667
#> 45 59_4_16 0.1666667
#> 46  59_8_3 0.1666667
#> 47  65_9_9 0.1666667
#> 
#> $occupancy_abundance
#>        otu   otu_occ     otu_rel fill
#> 1  59_8_22 1.0000000 0.121510673   no
#> 2  59_5_13 0.3333333 0.001642036   no
#> 3  59_8_12 0.3333333 0.018062397 core
#> 4  65_3_22 0.3333333 0.001642036   no
#> 5   65_5_1 0.3333333 0.001642036   no
#> 6  65_1_10 0.3333333 0.001642036   no
#> 7  65_7_12 1.0000000 0.027914614 core
#> 8   59_6_1 0.3333333 0.001642036   no
#> 9  65_2_17 1.0000000 0.019704433 core
#> 10 65_9_26 0.3333333 0.008210181   no
#> 11 65_5_18 0.3333333 0.001642036   no
#> 12 65_8_12 0.3333333 0.001642036   no
#> 13  65_9_1 0.3333333 0.003284072   no
#> 14 59_9_26 0.3333333 0.006568144   no
#> 15  9_6_28 0.3333333 0.001642036   no
#> 16  65_7_4 1.0000000 0.014778325 core
#> 17 65_4_26 0.3333333 0.001642036   no
#> 18 65_1_17 0.3333333 0.004926108   no
#> 19   9_4_3 0.3333333 0.001642036   no
#> 20  65_8_7 0.3333333 0.003284072   no
#> 21  59_4_5 0.3333333 0.001642036   no
#> 22  59_5_2 1.0000000 0.021346470 core
#> 23  65_6_2 0.6666667 0.018062397 core
#> 24 65_4_10 0.3333333 0.004926108   no
#> 25  9_7_25 0.6666667 0.006568144   no
#> 26   9_4_5 0.3333333 0.006568144   no
#> 27  65_1_8 0.3333333 0.001642036   no
#> 28  65_2_5 1.0000000 0.009852217 core
#> 29  59_7_6 1.0000000 0.151067323 core
#> 30 59_9_31 1.0000000 0.026272578 core
#> 31  9_4_13 0.3333333 0.001642036   no
#> 32   9_6_3 0.3333333 0.003284072   no
#> 33 59_9_18 0.3333333 0.004926108   no
#> 34 65_9_13 0.6666667 0.008210181   no
#> 35 65_3_18 0.6666667 0.013136289 core
#> 36 65_8_25 0.3333333 0.004926108   no
#> 37 65_8_29 0.3333333 0.001642036   no
#> 38 65_4_20 0.3333333 0.004926108   no
#> 39 59_9_17 0.6666667 0.013136289 core
#> 40   9_4_6 0.6666667 0.008210181   no
#> 41 59_4_25 0.3333333 0.011494253 core
#> 42 59_3_21 0.3333333 0.013136289 core
#> 43 59_4_16 0.3333333 0.001642036   no
#> 44  59_8_3 0.3333333 0.004926108   no
#> 45 59_5_19 1.0000000 0.073891626 core
#> 46  65_9_9 0.3333333 0.001642036   no
#> 47  59_2_6 1.0000000 0.336617406 core
#> 
#> $otu_table
#>          B  C   D
#> 59_8_22 55 16   3
#> 59_5_13  0  1   0
#> 59_8_12  0 11   0
#> 65_3_22  0  1   0
#> 65_5_1   0  0   1
#> 65_1_10  0  0   1
#> 65_7_12  5  2  10
#> 59_6_1   0  1   0
#> 65_2_17  7  1   4
#> 65_9_26  0  0   5
#> 65_5_18  0  0   1
#> 65_8_12  0  0   1
#> 65_9_1   0  0   2
#> 59_9_26  0  4   0
#> 9_6_28   1  0   0
#> 65_7_4   1  6   2
#> 65_4_26  0  0   1
#> 65_1_17  0  3   0
#> 9_4_3    1  0   0
#> 65_8_7   2  0   0
#> 59_4_5   0  1   0
#> 59_5_2   2  9   2
#> 65_6_2   0  7   4
#> 65_4_10  0  3   0
#> 9_7_25   1  0   3
#> 9_4_5    4  0   0
#> 65_1_8   0  1   0
#> 65_2_5   1  4   1
#> 59_7_6  39 36  17
#> 59_9_31  7  8   1
#> 9_4_13   1  0   0
#> 9_6_3    2  0   0
#> 59_9_18  3  0   0
#> 65_9_13  0  1   4
#> 65_3_18  2  6   0
#> 65_8_25  0  0   3
#> 65_8_29  0  0   1
#> 65_4_20  0  0   3
#> 59_9_17  1  0   7
#> 9_4_6    0  3   2
#> 59_4_25  0  7   0
#> 59_3_21  0  8   0
#> 59_4_16  0  1   0
#> 59_8_3   0  3   0
#> 59_5_19 15 25   5
#> 65_9_9   0  0   1
#> 59_2_6  53 34 118
#> 
#> $sample_metadata
#>   Sample Group SampleID
#> B      B     A        B
#> C      C     B        C
#> C      C     A        C
#> 
#> $taxonomy_table
#>          Kingdom         Phylum               Class              Order
#> 59_8_22  Archaea  Euryarchaeota         Bacteroidia Methanobacteriales
#> 59_5_13 Bacteria Proteobacteria Gammaproteobacteria      Bacteroidales
#> 59_8_12  Archaea  Euryarchaeota     Methanobacteria      Clostridiales
#> 65_3_22 Bacteria     Firmicutes Gammaproteobacteria Methanobacteriales
#> 65_5_1   Archaea     Firmicutes      Actinobacteria      Clostridiales
#> 65_1_10  Archaea Actinobacteria          Clostridia   Enterobacterales
#> 65_7_12 Bacteria  Bacteroidetes     Methanobacteria  Bifidobacteriales
#> 59_6_1   Archaea     Firmicutes     Methanobacteria      Bacteroidales
#> 65_2_17 Bacteria Actinobacteria         Bacteroidia  Bifidobacteriales
#> 65_9_26  Archaea     Firmicutes      Actinobacteria Methanobacteriales
#> 65_5_18  Archaea Actinobacteria          Clostridia      Bacteroidales
#> 65_8_12 Bacteria     Firmicutes     Methanobacteria   Enterobacterales
#> 65_9_1   Archaea  Bacteroidetes         Bacteroidia Methanobacteriales
#> 59_9_26 Bacteria Actinobacteria Gammaproteobacteria Methanobacteriales
#> 9_6_28   Archaea Proteobacteria          Clostridia      Bacteroidales
#> 65_7_4   Archaea Proteobacteria Gammaproteobacteria      Clostridiales
#> 65_4_26  Archaea  Bacteroidetes     Methanobacteria Methanobacteriales
#> 65_1_17  Archaea  Euryarchaeota Gammaproteobacteria   Enterobacterales
#> 9_4_3    Archaea  Euryarchaeota     Methanobacteria      Bacteroidales
#> 65_8_7  Bacteria     Firmicutes         Bacteroidia      Bacteroidales
#> 59_4_5   Archaea Actinobacteria         Bacteroidia      Bacteroidales
#> 59_5_2   Archaea Actinobacteria          Clostridia   Enterobacterales
#> 65_6_2   Archaea  Bacteroidetes      Actinobacteria      Bacteroidales
#> 65_4_10  Archaea Proteobacteria          Clostridia Methanobacteriales
#> 9_7_25  Bacteria  Bacteroidetes      Actinobacteria  Bifidobacteriales
#> 9_4_5   Bacteria  Bacteroidetes     Methanobacteria      Clostridiales
#> 65_1_8  Bacteria Actinobacteria      Actinobacteria      Clostridiales
#> 65_2_5   Archaea Actinobacteria          Clostridia  Bifidobacteriales
#> 59_7_6  Bacteria Proteobacteria Gammaproteobacteria   Enterobacterales
#> 59_9_31  Archaea  Euryarchaeota Gammaproteobacteria   Enterobacterales
#> 9_4_13  Bacteria Actinobacteria      Actinobacteria Methanobacteriales
#> 9_6_3    Archaea  Bacteroidetes          Clostridia  Bifidobacteriales
#> 59_9_18  Archaea Proteobacteria      Actinobacteria Methanobacteriales
#> 65_9_13 Bacteria  Bacteroidetes Gammaproteobacteria      Bacteroidales
#> 65_3_18  Archaea Proteobacteria          Clostridia Methanobacteriales
#> 65_8_25  Archaea     Firmicutes      Actinobacteria Methanobacteriales
#> 65_8_29 Bacteria  Euryarchaeota Gammaproteobacteria      Clostridiales
#> 65_4_20 Bacteria Actinobacteria     Methanobacteria   Enterobacterales
#> 59_9_17 Bacteria     Firmicutes     Methanobacteria      Bacteroidales
#> 9_4_6    Archaea  Bacteroidetes     Methanobacteria  Bifidobacteriales
#> 59_4_25 Bacteria     Firmicutes Gammaproteobacteria      Clostridiales
#> 59_3_21 Bacteria Proteobacteria Gammaproteobacteria      Clostridiales
#> 59_4_16 Bacteria  Bacteroidetes      Actinobacteria  Bifidobacteriales
#> 59_8_3  Bacteria     Firmicutes          Clostridia      Clostridiales
#> 59_5_19  Archaea     Firmicutes     Methanobacteria Methanobacteriales
#> 65_9_9   Archaea  Euryarchaeota      Actinobacteria      Bacteroidales
#> 59_2_6  Bacteria     Firmicutes         Bacteroidia      Clostridiales
#>                      Family              Genus                    Species
#> 59_8_22      Bacteroidaceae    Bifidobacterium           Escherichia coli
#> 59_5_13  Bifidobacteriaceae        Bacteroides           Blautia producta
#> 59_8_12     Lachnospiraceae        Escherichia           Blautia producta
#> 65_3_22 Methanobacteriaceae            Blautia           Escherichia coli
#> 65_5_1  Methanobacteriaceae Methanobrevibacter     Bifidobacterium longum
#> 65_1_10  Bifidobacteriaceae Methanobrevibacter           Blautia producta
#> 65_7_12  Enterobacteriaceae        Escherichia Methanobrevibacter smithii
#> 59_6_1   Bifidobacteriaceae Methanobrevibacter           Blautia producta
#> 65_2_17  Enterobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 65_9_26  Enterobacteriaceae        Bacteroides       Bacteroides fragilis
#> 65_5_18  Enterobacteriaceae        Bacteroides           Blautia producta
#> 65_8_12      Bacteroidaceae            Blautia Methanobrevibacter smithii
#> 65_9_1  Methanobacteriaceae Methanobrevibacter       Bacteroides fragilis
#> 59_9_26 Methanobacteriaceae            Blautia     Bifidobacterium longum
#> 9_6_28   Bifidobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 65_7_4  Methanobacteriaceae            Blautia           Escherichia coli
#> 65_4_26     Lachnospiraceae        Escherichia           Escherichia coli
#> 65_1_17      Bacteroidaceae            Blautia     Bifidobacterium longum
#> 9_4_3    Bifidobacteriaceae            Blautia           Escherichia coli
#> 65_8_7       Bacteroidaceae            Blautia Methanobrevibacter smithii
#> 59_4_5   Enterobacteriaceae            Blautia     Bifidobacterium longum
#> 59_5_2       Bacteroidaceae        Escherichia           Escherichia coli
#> 65_6_2   Enterobacteriaceae Methanobrevibacter           Blautia producta
#> 65_4_10  Enterobacteriaceae            Blautia           Escherichia coli
#> 9_7_25   Enterobacteriaceae Methanobrevibacter     Bifidobacterium longum
#> 9_4_5   Methanobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 65_1_8       Bacteroidaceae            Blautia       Bacteroides fragilis
#> 65_2_5       Bacteroidaceae Methanobrevibacter           Escherichia coli
#> 59_7_6   Enterobacteriaceae        Escherichia     Bifidobacterium longum
#> 59_9_31     Lachnospiraceae        Escherichia           Blautia producta
#> 9_4_13      Lachnospiraceae        Bacteroides           Blautia producta
#> 9_6_3    Bifidobacteriaceae            Blautia           Escherichia coli
#> 59_9_18  Bifidobacteriaceae        Bacteroides     Bifidobacterium longum
#> 65_9_13      Bacteroidaceae            Blautia           Escherichia coli
#> 65_3_18  Enterobacteriaceae        Escherichia Methanobrevibacter smithii
#> 65_8_25  Bifidobacteriaceae        Bacteroides       Bacteroides fragilis
#> 65_8_29  Enterobacteriaceae        Escherichia Methanobrevibacter smithii
#> 65_4_20 Methanobacteriaceae        Bacteroides           Escherichia coli
#> 59_9_17  Enterobacteriaceae        Escherichia Methanobrevibacter smithii
#> 9_4_6        Bacteroidaceae        Bacteroides           Blautia producta
#> 59_4_25  Bifidobacteriaceae        Escherichia Methanobrevibacter smithii
#> 59_3_21 Methanobacteriaceae        Bacteroides     Bifidobacterium longum
#> 59_4_16  Enterobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 59_8_3      Lachnospiraceae        Bacteroides           Escherichia coli
#> 59_5_19 Methanobacteriaceae        Escherichia       Bacteroides fragilis
#> 65_9_9      Lachnospiraceae    Bifidobacterium           Blautia producta
#> 59_2_6   Enterobacteriaceae            Blautia           Escherichia coli
```

## Output

The `extract_core()` function returns a list containing the following
elements:

- **core_otus**: A vector of core OTUs identified based on the specified
  threshold.

- **bray_curtis_ranked**: A data frame of Bray-Curtis dissimilarity
  rankings.

- **otu_rankings**: A data frame of OTU rankings based on occupancy and
  abundance.

- **occupancy_abundance**: A data frame of occupancy and abundance
  values for each OTU.

- **otu_table**: The OTU table used for analysis.

- **sample_metadata**: The sample metadata used for analysis.

- **taxonomy_table**: The taxonomy table used for analysis.
