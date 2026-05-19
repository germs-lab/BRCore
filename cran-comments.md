## v2.0.6
### R CMD check results 

0 errors | 0 warnings | 3 notes

### Fixes

* Relaxed version requirements for `phyloseq` and `stats4` to improve compatibility and avoid 
CRAN error in `r-oldrel-macos-x86_64`.

### Fixes

* Relaxed version requirements for `phyloseq` and `stats4` to improve compatibility and avoid 
CRAN error in `r-oldrel-macos-x86_64`.

### Notes

1. **CRAN incoming feasibility**:
   - `Running aspell failed`: The `aspell` English dictionary is unavailable in
     the local check environment (`/usr/lib/aspell-0.60/en_US`). This is a
     local environment issue and does not affect the package.

2. **Examples with CPU (user + system) or elapsed time > 5s**
                   user system elapsed
  update_otu_table 5.01  0.034   5.055
2. **Examples with CPU (user + system) or elapsed time > 5s**
                   user system elapsed
  update_otu_table 5.01  0.034   5.055

3. **HTML version of manual**: `V8` package unavailable in the local check
   environment, skipping math rendering check. This is a local environment
   issue.




