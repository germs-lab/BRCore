# BRCore 2.0.0
Date: 2026-04-24

## Breaking Changes

* [de1fa42](https://github.com/germs-lab/BRCore/commit/de1fa425d7466f7e652136b949884f795d878e0c) `multi_rarefy()`: `.as_array` logical argument replaced with `.as` character argument (e.g. `"list"`, `"array"`). Single-iteration handling improved to support all `.as` formats and seeds.
* [fbb3b90](https://github.com/germs-lab/BRCore/commit/fbb3b907b7010a3c64e6efb096b5d797b8f90fd9) `identify_core()`: output column renamed from `Index` to `rank`.
* [fa7387e](https://github.com/germs-lab/BRCore/commit/fa7387e6e74e9cb4ac5dafe4f993288717a8ba54) `find_core()` removed; ground truth logic moved to `debugging` branch.
* [28b8821](https://github.com/germs-lab/BRCore/commit/28b8821a49e438c0a0bd5d14dd878689b9b97bff) and [4f5a9c9](https://github.com/germs-lab/BRCore/commit/4f5a9c91ea4e5972d0937680a63704a3d2a0c9ab) `identify_core()` major rewrite: `vegan::avgdist()` replaced with `vegan::vegdist()` via new `.mean_bc_over_iters()` internal helper; `rarefied_list` is now optional and requires ≥ 2 iterations when provided. This fixes incorrect results on unrarefied datasets caused by how `vegan::avgdist()` handles zeros in a matrix of ones.
* [9f3b3ae](https://github.com/germs-lab/BRCore/commit/9f3b3ae29b999b4691c4e0fffb9b9e358063e389) `plot_identified_core()` now returns a named list with `$df` and `$plot` instead of a bare plot object.
* [3a590aa](https://github.com/germs-lab/BRCore/commit/3a590aa0b3a1ebccf868ae1186804d64ccab7f21) General update/refactor of tests to include new function parameters and logic. Regenerated `tests/testthat/test_sets/test_vignette_data.rda` reference dataset.
* [2677cd2](https://github.com/germs-lab/BRCore/commit/2677cd23ecf5a35513fd6a40b5d6ee4a3bc76919) Fixed binomial neutral model fitting: `pbinom()` now receives `N.int`, an integer `size` argument (`round(N)`) instead of a floating-point mean. This changes model output relative to prior versions. We also added comprehensive messaging to inform user of species (`N` or `N.int`) used in model fitting.
* [200100b](https://github.com/germs-lab/BRCore/commit/200100b8306371cb1ca7aacf130f97bf056f66fb) `multi_rarefy()` internal logic: Replaced `.single_rarefy()` with internal `vegan::rrarefy()` engine.
* [34c41ca](https://github.com/germs-lab/BRCore/commit/34c41caa1477aa29dfa8dfde02c7d6324082c96c) and [114468c](https://github.com/germs-lab/BRCore/commit/114468cf60b6276e5cd9a3829ae43874a8e087bf) Refactored `multi_rarefy()` to output 3D array or list. This handling increased speed of computation and eliminated the need for parallelization. CLI messages were improved to handle these new types and present summary statistics and information to users.

## New Features

* [dbdf9fe](https://github.com/germs-lab/BRCore/commit/dbdf9fecd5dbecde83bba3a72038ffcbfa267e71) Added `.brcore_theme()` internal helper (`brcore_theme.R`) to unify plot styling (borders, title sizes, viridis palettes) across all plotting functions.
* `plot_identified_core()` gains an optional `dataset_name` parameter for plot titles.

## Bug Fixes

* [3225bd4](https://github.com/germs-lab/BRCore/commit/3225bd4c28180ded6e36bc0e7fbdfa2a2496f7cd) `identify_core()`: fixed `proportionBC` normalisation to use `max()` instead of `last()`.
* [b3824ab](https://github.com/germs-lab/BRCore/commit/b3824ab272ec16b65c0136e10e65e5925ff74474) `identify_core()`: fixed BC ranking and pair alignment using unique time points.
* [7439ed3](https://github.com/germs-lab/BRCore/commit/7439ed31dc496d60165f7a250e957b6eabfc7f82) `plot_identified_core()`: fixed deprecated `size` → `linewidth` in `panel.border`.
* [5c2c57e](https://github.com/germs-lab/BRCore/commit/5c2c57eaa574b61bfa8045ef4614d64413680a68) and [bfd74da](https://github.com/germs-lab/BRCore/commit/bfd74da34cd63c1ef59c491679fda49299cd7701) Added `plot_variance_propagation()` to plot results from `multi_rarefy()`.
* [79153ca](https://github.com/germs-lab/BRCore/commit/79153cabcdd93f44080566ea09d18f194afbae9d) `update_otu_table()` has new parameter `iteration` to handle results from `multi_rarefy()`.
* [e4ac648](https://github.com/germs-lab/BRCore/commit/e4ac648420ceaa6bb3d5d946e8f7170029e6f615) `identify_core()` now has a progress bar powered by `cli::cli_progress_bar()` in `.calculate_bc()`.

## Smaller Changes

* [36285d2](https://github.com/germs-lab/BRCore/commit/36285d232d595702ce9846566b40957febd4777a) Tidy eval compliance enforced throughout using `.data$` pronouns and quoted column names; `globals.R` significantly trimmed.
* Vignette and README updated to reflect the required multi-iteration workflow for `identify_core()`.
* Test suite updated for `identify_core`, `multi_rarefy`, `plot_identified_core`, `plot_variance_propagation`, `update_otu_table`, and the full vignette workflow. Test data regenerated (`test_vignette_data.rda` grew ~250 KB).

* Vignette updated to include `multi_rarefy()` and `plot_variance_propagation()` as part of the Data Exploration and Parameter Selection section. This new workflow allows users to determine the adequate sequence read depth for their data and input that into `identify_core()`.

**See PR [#97](https://github.com/germs-lab/BRCore/pull/97) and [#95](https://github.com/germs-lab/BRCore/pull/95) for more details**

## Contributors
@jibarozzo @Gian77

# BRCore 1.0.2
Date: 2026-02-25

## Bug Fixes

* Fixed `multi_rarefy()` reproducibility across platforms using deterministic iteration seeds for parallel processing and adding floating-point tolerance in rarefaction.
    - **Commits:** [`0422317`](https://github.com/germs-lab/BRCore/commit/042231709ec858a2f5869886aba17f82bef2c938), [`2341154`](https://github.com/germs-lab/BRCore/commit/23411546814a8bdf9c1346f22e42f1820b77c371), [`a42d7a9`](https://github.com/germs-lab/BRCore/commit/a42d7a96e1931d99825fce1c66a4383d51e7b5c6)

* Added floating point tolerance in `identify_core()` rarefaction validation
    - **Commits:** [`a687fa3`](https://github.com/germs-lab/BRCore/commit/a687fa38327dab62b9c7f666d25a59dc75722c05)

## Improvements

* Changed rarefaction engine `vegan::rrarefy()` to custom `.single_rarefy()` reducing reliance on external dependencies. 
    - **Commits:** [`d366614`](https://github.com/germs-lab/BRCore/commit/d366614b0a715ff3ddf2cdfed18ebe93d872962b), [`76fbec2`](https://github.com/germs-lab/BRCore/commit/76fbec2b6e3e4f72e98ae6fc5bd2e7908cfd1045)
* Enhanced CLI output in `multi_rarefy()`: improved input validation for `multi_rarefy()` parameters
    - **Commits::** [`921dc51`](https://github.com/germs-lab/BRCore/commit/921dc51496fe22e65ac4a56eff11b53c515dc259)
* Added comprehensive end-to-end workflow tests that replicate the vignette workflow. 

    - **Notable Commits:** [`3080641`](https://github.com/germs-lab/BRCore/commit/308064172ebb1b6c2f7a82374c78ba4b4b976bfb), [`622b2a7`](https://github.com/germs-lab/BRCore/commit/622b2a7bac3673cd42216935f91ce3d5133f7d92), [`3877398`](https://github.com/germs-lab/BRCore/commit/3877398ad5817f4c7faaffa49bacde4d779fca93)

## Testing

* Added 4 comprehensive test suites validating complete analysis pipeline
* Reference data regeneration workflow documented
    - **Notable Commits:** [`924ca44`](https://github.com/germs-lab/BRCore/commit/924ca44c4ee3b7bfce35e6a3f05ae3b59d5e069c), [`cc8d81b`](https://github.com/germs-lab/BRCore/commit/cc8d81b32bdc01165374d42630244aea4b4116bd), [`2341154`](https://github.com/germs-lab/BRCore/commit/23411546814a8bdf9c1346f22e42f1820b77c371) 
* All tests use appropriate floating-point tolerance

**See PR #84 for more details**

## Contributors
@Gian77 @jibarozzo @Copilot
