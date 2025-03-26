# Test suite
test_that("extract_core() works with esophagus_with_tax dataset", {
  # Load the esophagus dataset
  data(esophagus, package = "phyloseq")
  
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
    Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE),
    # Random groups
    row.names = sample_names(esophagus)
  )
  
  # Add the taxonomy table to the esophagus dataset
  esophagus_with_tax <- merge_phyloseq(esophagus, tax_table, sample_data(sample_data))
  
  
  # Load expected
  load(here::here("tests/testthat/expected_extract_core.rda"))
  
  # Run the function
  test_core <- extract_core(
    esophagus_with_tax,
    Var = "site",
    method = "increase",
    increase_value = 2
  )
  testthat::expect_identical(test_core, expected_extract_core)
  
  # Test 1: Check output structure
  expect_named(
    test_core,
    c(
      "core_otus",
      "bray_curtis_ranked",
      "otu_rankings",
      "occupancy_abundance",
      "otu_table",
      "sample_metadata",
      "taxonomy_table"
    )
  )
  
  # Test 2: Check core_otus
  expect_true(is.character(test_core$core_otus))
  expect_true(length(test_core$core_otus) > 0)
  
  # Test 3: Check bray_curtis_ranked
  expect_s3_class(test_core$bray_curtis_ranked, "data.frame")
  expect_true(all(
    c("rank", "MeanBC", "proportionBC", "IncreaseBC") %in% names(test_core$bray_curtis_ranked)
  ))
  
  # Test 4: Check otu_rankings
  expect_s3_class(test_core$otu_rankings, "data.frame")
  expect_true(all(c("otu", "rank") %in% names(test_core$otu_rankings)))
  
  # Test 5: Check occupancy_abundance
  expect_s3_class(test_core$occupancy_abundance, "data.frame")
  expect_true(all(
    c("otu", "otu_occ", "otu_rel") %in% names(test_core$occupancy_abundance)
  ))
  
  # Test 6: Check otu_table
  expect_true(is.matrix(test_core$otu_table))
  expect_equal(nrow(test_core$otu_table),
               nrow(expected_extract_core$otu_table))
  
  # Test 7: Check sample_metadata
  expect_s3_class(test_core$sample_metadata, "data.frame")
  expect_true(all(
    c("Sample", "SampleID") %in% names(test_core$sample_metadata)
  ))
  
  # Test 8: Check taxonomy_table
  expect_s3_class(test_core$taxonomy_table, "data.frame")
  expect_true(all(
    c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    ) %in% names(test_core$taxonomy_table)
  ))
  
  # Test 9: Check for NA/NaN in bray_curtis_ranked
  expect_false(any(is.na(test_core$bray_curtis_ranked)))
  
  # Test 10: Check for NA/NaN in otu_rankings
  expect_false(any(is.na(test_core$otu_rankings)))
})
