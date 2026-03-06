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
