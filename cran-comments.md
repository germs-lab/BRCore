## R CMD check results (v2.0.1 resubmission)

0 errors | 0 warnings | 2 notes

### Notes

1. **CRAN incoming feasibility**:
   - `New submission`: Expected for new submission or re-submission.
   - `Running aspell failed`: The `aspell` English dictionary is unavailable in
     the local check environment (`/usr/lib/aspell-0.60/en_US`). This is a
     local environment issue and does not affect the package.

2. **HTML version of manual**: `V8` package unavailable in the local check
   environment, skipping math rendering check. This is a local environment
   issue.

### Resubmission notes (v2.0.1)

This is a resubmission of v2.0.0, which was rejected by the CRAN auto-check
service (win-builder, r-devel-windows-x86_64). The following issues were
identified and resolved:

1. **URL in DESCRIPTION**: Changed `http://` to `https://` for the package
   website URL.

2. **Description field**: Reworded to not start with "This package" per CRAN
   policy.

3. **Possibly misspelled words** (`BRCore`, `Bioenergy`, `microbiome`, and `microbiomes`): These are intentional. `BRCore` is the package name and `Bioenergy` , `microbiome`, and `microbiomes` are domain-specific terms.

4. **Invalid file URIs in README.md**: Relative links to `CONTRIBUTING.md` and
   `CODE_OF_CONDUCT.md` replaced with full GitHub URLs.

### win-builder notes

The win-builder (r-devel-windows-x86_64) reported errors for missing packages
(`stringi`, `jsonlite`, `xfun`, `Formula`). These are indirect dependencies of
declared imports (`phyloseq`, `Hmisc`, `knitr`) that were absent in the
win-builder guest environment. The CRAN incoming checks on both Windows and
Debian passed with only NOTEs, confirming these are not direct imports in
package code.

* This is a new submission to CRAN.
