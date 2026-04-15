# tests/testthat/setup.R

rda_path <- test_path("test_sets/test_vignette_data.rda")

# Only regenerate if the .rda is missing or if running on CI
if (!file.exists(rda_path)) {
  source(test_path("test_sets/regenerate_test_data.R"))
}

if (nzchar(Sys.getenv("CI"))) {
  rda_path <- file.path(tempdir(), "test_sets/test_vignette_data.rda")
  options(brcore_test_data_path = rda_path)
  source(test_path("test_sets/regenerate_test_data.R"), local = TRUE)
}
