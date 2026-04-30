## Resubmission notes (v2.0.3)

This is a resubmission. Previous submissions:

- v2.0.2: Passed auto-checks; requested resubmission by manual inspection
- v2.0.1: Passed auto-checks; retracted by maintainer due to typos
- v2.0.0: Rejected by CRAN auto-check service

**The following issues from v2.0.2 were resolved:**

1. **DOI** links in the documentation were updated to use the correct format (`<doi:...>`).
2. Minor typos and formatting issues in the documentation were corrected.
3. Examples were unwrapped from `\dontrun` and `\donttest`. Only example for `identify_core()` remains wrapped due to long execution time (> 5s).
4. We removed `options(warn=-1)` from `scnm.fit()` as it was not needed and remained from a debugging phase.
5. An appropriate example for `scnm.fit()` was included.



### R CMD check results 

0 errors | 0 warnings | 3 notes

### Notes

1. **CRAN incoming feasibility**:
   - `New submission`: Expected for new submission or re-submission.
   - `Running aspell failed`: The `aspell` English dictionary is unavailable in
     the local check environment (`/usr/lib/aspell-0.60/en_US`). This is a
     local environment issue and does not affect the package.

2. **Checking for future file timestamps**: Unable to verify current time. This is a
     local environment issue and does not affect the package.

3. **HTML version of manual**: `V8` package unavailable in the local check
   environment, skipping math rendering check. This is a local environment
   issue.

---

### win-builder notes

The win-builder (r-devel-windows-x86_64) reported errors for missing packages
(`stringi`, `jsonlite`, `xfun`, `Formula`). These are indirect dependencies of
declared imports (`phyloseq`, `Hmisc`, `knitr`) that were absent in the
win-builder guest environment. The CRAN incoming checks on both Windows and
Debian passed with only NOTEs, confirming these are not direct imports in
package code.
 


