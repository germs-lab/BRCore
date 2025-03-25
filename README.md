
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
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#> 
#> → Will update 1 package.
#> → The package (0 B) is cached.
#> + BRCore 0.0.0.9000 → 0.0.0.9000 [bld][cmp] (GitHub: 51d757d)
#> ✔ All system requirements are already installed.
#> 
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got BRCore 0.0.0.9000 (source) (740.97 kB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [142 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcurl4-openssl-dev libssl-dev make libglpk-dev libxml2-dev pandoc libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.5).
#> make is already the newest version (4.3-4.1build2).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.2).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 54 not upgraded.
#> ℹ Packaging BRCore 0.0.0.9000
#> ✔ Packaged BRCore 0.0.0.9000 (564ms)
#> ℹ Building BRCore 0.0.0.9000
#> ✔ Built BRCore 0.0.0.9000 (8.8s)
#> ✔ Installed BRCore 0.0.0.9000 (github::germs-lab/BRCore@51d757d) (1s)
#> ✔ 1 pkg + 114 deps: kept 111, upd 1, dld 1 (NA B) [17.9s]
```

## Step 1: Load Required Libraries and Data

``` r
library(phyloseq)
library(BRCore)

# Load the esophagus dataset
data(esophagus, package = "phyloseq")
```

## Step 2: Create Mock Taxonomy Table

``` r
# Get the taxa names from the esophagus dataset
taxa_names <- taxa_names(esophagus)

# Define realistic taxonomy levels
kingdoms <- c("Bacteria", "Archaea")
phyla <- c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria", "Euryarchaeota")
classes <- c("Clostridia", "Bacteroidia", "Gammaproteobacteria", "Actinobacteria", "Methanobacteria")
orders <- c("Clostridiales", "Bacteroidales", "Enterobacterales", "Bifidobacteriales", "Methanobacteriales")
families <- c("Lachnospiraceae", "Bacteroidaceae", "Enterobacteriaceae", "Bifidobacteriaceae", "Methanobacteriaceae")
genera <- c("Blautia", "Bacteroides", "Escherichia", "Bifidobacterium", "Methanobrevibacter")
species <- c("Blautia producta", "Bacteroides fragilis", "Escherichia coli", "Bifidobacterium longum", "Methanobrevibacter smithii")

# Create a mock taxonomy table
mock_taxonomy_table <- matrix(
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
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
)

# Convert to a taxonomy table object
tax_table <- tax_table(mock_taxonomy_table)
```

## Step 3: Add Sample Metadata

``` r
# Add sample metadata to the esophagus dataset
sample_data <- data.frame(
  Sample = sample_names(esophagus),
  Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE), # Random groups
  row.names = sample_names(esophagus)
)

# Add the taxonomy table and sample metadata to the esophagus dataset
esophagus_with_tax <- merge_phyloseq(esophagus, tax_table, sample_data(sample_data))
```

## Step 4: Run extract_core()

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
#> ℹ Calculating BC dissimilarity based on ranked OTUs, starting at 2025-03-21 15:54:35.26651
#> ✔ BC ranks done!
#> Joining with `by = join_by(rank)`
#> Joining with `by = join_by(rank, IncreaseBC)`
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
#> C      C     A        C
#> D      D     B        D
#> 
#> $taxonomy_table
#>          Kingdom         Phylum               Class              Order
#> 59_8_22  Archaea  Euryarchaeota     Methanobacteria  Bifidobacteriales
#> 59_5_13 Bacteria  Euryarchaeota Gammaproteobacteria   Enterobacterales
#> 59_8_12 Bacteria     Firmicutes          Clostridia   Enterobacterales
#> 65_3_22 Bacteria  Bacteroidetes         Bacteroidia      Bacteroidales
#> 65_5_1   Archaea     Firmicutes         Bacteroidia      Clostridiales
#> 65_1_10  Archaea  Euryarchaeota     Methanobacteria      Clostridiales
#> 65_7_12  Archaea Proteobacteria          Clostridia      Bacteroidales
#> 59_6_1  Bacteria  Bacteroidetes          Clostridia Methanobacteriales
#> 65_2_17  Archaea  Bacteroidetes      Actinobacteria  Bifidobacteriales
#> 65_9_26 Bacteria  Euryarchaeota     Methanobacteria  Bifidobacteriales
#> 65_5_18 Bacteria Proteobacteria         Bacteroidia   Enterobacterales
#> 65_8_12 Bacteria     Firmicutes Gammaproteobacteria  Bifidobacteriales
#> 65_9_1  Bacteria  Bacteroidetes      Actinobacteria  Bifidobacteriales
#> 59_9_26  Archaea  Bacteroidetes          Clostridia Methanobacteriales
#> 9_6_28  Bacteria Actinobacteria     Methanobacteria  Bifidobacteriales
#> 65_7_4  Bacteria     Firmicutes     Methanobacteria   Enterobacterales
#> 65_4_26  Archaea Actinobacteria         Bacteroidia Methanobacteriales
#> 65_1_17  Archaea Proteobacteria     Methanobacteria   Enterobacterales
#> 9_4_3   Bacteria     Firmicutes Gammaproteobacteria      Clostridiales
#> 65_8_7   Archaea  Euryarchaeota      Actinobacteria      Clostridiales
#> 59_4_5   Archaea Proteobacteria          Clostridia Methanobacteriales
#> 59_5_2  Bacteria  Euryarchaeota     Methanobacteria Methanobacteriales
#> 65_6_2   Archaea  Euryarchaeota Gammaproteobacteria   Enterobacterales
#> 65_4_10  Archaea Actinobacteria      Actinobacteria      Bacteroidales
#> 9_7_25  Bacteria     Firmicutes      Actinobacteria Methanobacteriales
#> 9_4_5   Bacteria     Firmicutes      Actinobacteria   Enterobacterales
#> 65_1_8  Bacteria Proteobacteria Gammaproteobacteria      Clostridiales
#> 65_2_5  Bacteria     Firmicutes          Clostridia   Enterobacterales
#> 59_7_6  Bacteria  Bacteroidetes      Actinobacteria      Bacteroidales
#> 59_9_31  Archaea  Euryarchaeota      Actinobacteria   Enterobacterales
#> 9_4_13  Bacteria Proteobacteria         Bacteroidia      Clostridiales
#> 9_6_3    Archaea Actinobacteria         Bacteroidia      Clostridiales
#> 59_9_18 Bacteria  Bacteroidetes         Bacteroidia  Bifidobacteriales
#> 65_9_13  Archaea     Firmicutes         Bacteroidia      Clostridiales
#> 65_3_18 Bacteria     Firmicutes     Methanobacteria   Enterobacterales
#> 65_8_25  Archaea  Euryarchaeota Gammaproteobacteria  Bifidobacteriales
#> 65_8_29  Archaea  Bacteroidetes         Bacteroidia   Enterobacterales
#> 65_4_20  Archaea Actinobacteria     Methanobacteria  Bifidobacteriales
#> 59_9_17 Bacteria Actinobacteria     Methanobacteria  Bifidobacteriales
#> 9_4_6    Archaea  Bacteroidetes      Actinobacteria      Clostridiales
#> 59_4_25 Bacteria Actinobacteria          Clostridia   Enterobacterales
#> 59_3_21  Archaea  Bacteroidetes         Bacteroidia  Bifidobacteriales
#> 59_4_16 Bacteria     Firmicutes     Methanobacteria   Enterobacterales
#> 59_8_3   Archaea Proteobacteria Gammaproteobacteria Methanobacteriales
#> 59_5_19 Bacteria  Euryarchaeota     Methanobacteria  Bifidobacteriales
#> 65_9_9  Bacteria Proteobacteria     Methanobacteria      Clostridiales
#> 59_2_6   Archaea  Bacteroidetes     Methanobacteria  Bifidobacteriales
#>                      Family              Genus                    Species
#> 59_8_22 Methanobacteriaceae Methanobrevibacter           Escherichia coli
#> 59_5_13  Enterobacteriaceae        Escherichia Methanobrevibacter smithii
#> 59_8_12     Lachnospiraceae Methanobrevibacter           Escherichia coli
#> 65_3_22     Lachnospiraceae        Bacteroides       Bacteroides fragilis
#> 65_5_1  Methanobacteriaceae    Bifidobacterium           Blautia producta
#> 65_1_10     Lachnospiraceae        Escherichia       Bacteroides fragilis
#> 65_7_12  Enterobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 59_6_1  Methanobacteriaceae Methanobrevibacter       Bacteroides fragilis
#> 65_2_17     Lachnospiraceae        Escherichia           Escherichia coli
#> 65_9_26 Methanobacteriaceae        Escherichia           Blautia producta
#> 65_5_18  Bifidobacteriaceae            Blautia           Escherichia coli
#> 65_8_12     Lachnospiraceae            Blautia     Bifidobacterium longum
#> 65_9_1   Bifidobacteriaceae    Bifidobacterium       Bacteroides fragilis
#> 59_9_26  Bifidobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 9_6_28   Bifidobacteriaceae Methanobrevibacter           Blautia producta
#> 65_7_4   Enterobacteriaceae Methanobrevibacter       Bacteroides fragilis
#> 65_4_26  Enterobacteriaceae Methanobrevibacter Methanobrevibacter smithii
#> 65_1_17  Enterobacteriaceae        Bacteroides           Escherichia coli
#> 9_4_3        Bacteroidaceae Methanobrevibacter           Escherichia coli
#> 65_8_7       Bacteroidaceae        Bacteroides     Bifidobacterium longum
#> 59_4_5  Methanobacteriaceae        Bacteroides       Bacteroides fragilis
#> 59_5_2       Bacteroidaceae        Bacteroides Methanobrevibacter smithii
#> 65_6_2       Bacteroidaceae            Blautia           Escherichia coli
#> 65_4_10  Enterobacteriaceae    Bifidobacterium     Bifidobacterium longum
#> 9_7_25   Enterobacteriaceae            Blautia           Escherichia coli
#> 9_4_5        Bacteroidaceae        Bacteroides           Blautia producta
#> 65_1_8   Enterobacteriaceae        Escherichia           Blautia producta
#> 65_2_5   Enterobacteriaceae            Blautia           Blautia producta
#> 59_7_6  Methanobacteriaceae        Escherichia     Bifidobacterium longum
#> 59_9_31 Methanobacteriaceae            Blautia       Bacteroides fragilis
#> 9_4_13   Bifidobacteriaceae    Bifidobacterium           Blautia producta
#> 9_6_3       Lachnospiraceae        Bacteroides           Escherichia coli
#> 59_9_18  Enterobacteriaceae Methanobrevibacter     Bifidobacterium longum
#> 65_9_13 Methanobacteriaceae        Escherichia           Blautia producta
#> 65_3_18      Bacteroidaceae        Escherichia           Escherichia coli
#> 65_8_25  Enterobacteriaceae            Blautia           Blautia producta
#> 65_8_29     Lachnospiraceae            Blautia       Bacteroides fragilis
#> 65_4_20     Lachnospiraceae        Bacteroides       Bacteroides fragilis
#> 59_9_17     Lachnospiraceae        Bacteroides       Bacteroides fragilis
#> 9_4_6   Methanobacteriaceae            Blautia       Bacteroides fragilis
#> 59_4_25 Methanobacteriaceae            Blautia     Bifidobacterium longum
#> 59_3_21      Bacteroidaceae        Bacteroides Methanobrevibacter smithii
#> 59_4_16      Bacteroidaceae        Bacteroides       Bacteroides fragilis
#> 59_8_3  Methanobacteriaceae        Bacteroides     Bifidobacterium longum
#> 59_5_19     Lachnospiraceae        Escherichia Methanobrevibacter smithii
#> 65_9_9       Bacteroidaceae    Bifidobacterium     Bifidobacterium longum
#> 59_2_6   Bifidobacteriaceae            Blautia           Escherichia coli
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
