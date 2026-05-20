# tests/testthat/setup.R

rda_path <- test_path("test_sets/test_vignette_data.rda")

# Only regenerate if the .rda is missing or if running on CI
is_check <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
is_not_cran <- identical(Sys.getenv("NOT_CRAN"), "true")
is_ci <- nzchar(Sys.getenv("CI"))

needs <- c("devtools", "phyloseq", "tidyverse", "viridis")
can_regenerate <- all(vapply(
  needs,
  requireNamespace,
  logical(1),
  quietly = TRUE
))


# Does not regenerate on CRAN or during R CMD check, but does regenerate during devtools::test() if the .rda is missing
if (is_not_cran && !is_check) {
  if (!file.exists(rda_path)) {
    source(test_path("test_sets/regenerate_test_data.R"))
  }
}

# Regenerate test data if running on CI
if (is_ci && can_regenerate) {
  rda_path <- file.path(
    tempdir(),
    test_path("test_sets/test_vignette_data.rda")
  )
  options(brcore_test_data_path = rda_path)
  source(test_path("test_sets/regenerate_test_data.R"), local = TRUE)
}
