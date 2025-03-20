# Test suite
test_that("filter_core() works with esophagus dataset", {
  # Load the esophagus dataset
  data(esophagus, package = "phyloseq")
  
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
  
  # Add sample metadata to the esophagus dataset (required for calculate_nmds)
  sample_data <- data.frame(
    Sample = sample_names(esophagus),
    Group = sample(c("A", "B"), nsamples(esophagus), replace = TRUE), # Random groups
    row.names = sample_names(esophagus)
  )
  
  # Add the taxonomy table to the esophagus dataset
  esophagus_with_tax <- merge_phyloseq(esophagus, tax_table, sample_data(sample_data))
  
  
  # Test 1: Basic functionality
  results <- filter_core(esophagus_with_tax, threshold = 0.6, as = "rows")

  expect_named(results, c("physeq_high_occ", "physeq_low_occ"))
  expect_s4_class(results$physeq_high_occ, "phyloseq")
  expect_s4_class(results$physeq_low_occ, "phyloseq")

  # Test 2: Check high occurrence ASVs
  high_occ_asvs <- taxa_names(results$physeq_high_occ)
  expect_true(length(high_occ_asvs) > 0)
  expect_true(all(rowSums(otu_table(results$physeq_high_occ) > 0) >= 2)) # 60% of 3 samples

  # Test 3: Check low occurrence ASVs
  low_occ_asvs <- taxa_names(results$physeq_low_occ)
  expect_true(length(low_occ_asvs) > 0)
  expect_true(all(rowSums(otu_table(results$physeq_low_occ) > 0) < 2)) # Less than 60% of 10 samples

  # Test 4: Edge case - threshold = 0 (all ASVs are high occurrence)
  expect_error(filter_core(esophagus_with_tax, threshold = 0.0, as = "rows"))

  # Test 5: Edge case - threshold = 1 (all ASVs are low occurrence)
  results_threshold_1 <- filter_core(esophagus_with_tax, threshold = 1, as = "rows")
  expect_equal(ntaxa(results_threshold_1$physeq_high_occ), 14)
  expect_equal(ntaxa(results_threshold_1$physeq_low_occ), 44)

  # Test 6: Invalid threshold
  expect_error(filter_core(esophagus_with_tax, threshold = -1, as = "rows"), "must be between 0 and 1")
  expect_error(filter_core(esophagus_with_tax, threshold = 2, as = "rows"), "must be between 0 and 1")

  # Test 7: Invalid phyloseq object
  expect_error(filter_core("not_a_phyloseq_object", threshold = 0.6, as = "rows"), "must be a phyloseq object")

  # Test 8: Invalid 'as' argument
  expect_error(
    filter_core(esophagus_with_tax, threshold = 0.6, as = "invalid"),
    "Arguments 'as' must be a string:'rows', 'cols' or 'columns'!"
  )

  # Test 9: Check taxa orientation
  results_rows <- filter_core(esophagus_with_tax, threshold = 0.6, as = "rows")
  results_cols <- filter_core(esophagus_with_tax, threshold = 0.6, as = "columns")

  expect_true(taxa_are_rows(results_rows$physeq_high_occ))
  expect_true(taxa_are_rows(results_rows$physeq_low_occ))
  expect_false(taxa_are_rows(results_cols$physeq_high_occ))
  expect_false(taxa_are_rows(results_cols$physeq_low_occ))
})
