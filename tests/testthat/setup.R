# tests/testthat/setup.R

rda_path <- test_path("test_sets/test_vignette_data.rda")

# Only regenerate if the .rda is missing or if running on CI
if (!file.exists(rda_path) || nzchar(Sys.getenv("CI"))) {
  source(test_path("test_sets/regenerate_test_data.R"))
}
