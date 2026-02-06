# Cross-Platform Random Number Generation (RNG) Reproducibility

## Problem Statement

The `multi_rarefy()` function was producing different results on macOS vs Linux platforms even when using the same seed value. This issue manifested in the vignette workflow tests where rarefied OTU tables differed between operating systems.

## Root Causes

### 1. Default RNG Algorithm Differences
- R's default RNG algorithm can vary between platforms and R versions
- Without explicit specification, different operating systems may use different algorithms
- Similar issues documented in:
  - [ranger package issue #533](https://github.com/imbs-hl/ranger/issues/533)
  - [StackOverflow discussion](https://stackoverflow.com/questions/48626086/same-seed-different-os-different-random-numbers-in-r)

### 2. Parallel RNG Challenges
- Manual seed offsetting (`set.seed(set_seed + i)`) is not the recommended approach
- Different thread counts can produce different results
- No guarantee of independent streams across workers

### 3. Potential BLAS Library Differences
- Matrix operations (e.g., `as.matrix()`) might use different BLAS libraries
- Different BLAS implementations (OpenBLAS, MKR, Apple Accelerate) could theoretically affect floating-point operations
- However, for the rarefaction use case, this is less likely to be the primary issue

## Solution Implemented

### Primary Fix: Explicit RNG Algorithm

```r
# Set explicit RNG algorithm for cross-platform reproducibility
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(set_seed)
```

**Why Mersenne-Twister?**
- Most widely used and tested RNG algorithm in R
- Default on most platforms (ensuring we're explicit)
- Excellent statistical properties
- Cross-platform consistency

### Parallel RNG: clusterSetRNGStream

```r
# Set reproducible RNG streams for parallel workers
clusterSetRNGStream(cl, iseed = set_seed)
```

**Advantages:**
- Uses L'Ecuyer-CMRG algorithm designed for parallel computing
- Provides independent random number streams for each worker
- Results are reproducible regardless of number of threads
- Recommended by R's parallel package documentation
- Eliminates manual seed manipulation

### RNG State Management

```r
# Save current RNG state to restore later
old_rng_kind <- RNGkind()
on.exit(RNGkind(old_rng_kind[1], old_rng_kind[2], old_rng_kind[3]), add = TRUE)
```

**Why restore?**
- Avoids side effects on user's global RNG settings
- Good practice for functions that modify global state
- Uses `on.exit()` to ensure restoration even if function errors

## Alternative Solutions Considered

### 1. Force Single-Threaded Execution
- **Pros**: Eliminates parallel RNG complexity
- **Cons**: Significant performance degradation, not scalable
- **Status**: Already implemented in tests via CI detection, but not ideal for production

### 2. Different Parallel Backend
- **Pros**: Could use different parallelization strategy
- **Cons**: More complex, may not solve cross-platform issues
- **Status**: Not needed with proper RNG stream management

### 3. BLAS Standardization
- **Pros**: Could eliminate floating-point variation
- **Cons**: Not user-controllable, unlikely to be the main issue
- **Status**: Not the primary cause for discrete random sampling

## Testing Strategy

### Unit Tests
- Existing `test-multi_rarefy.R` validates basic functionality
- Single-threaded execution for deterministic results

### Vignette Workflow Tests
- `test-vignette-workflow.R` compares against reference data
- Tests run with `threads = 1` on CI for reproducibility
- Validates entire analysis pipeline

### Cross-Platform Validation
- CI runs on both macOS and Linux (ubuntu-latest)
- Reference data generated on one platform should match results on others

## Implementation Notes

### Function Changes
- Added RNG algorithm specification
- Implemented proper parallel RNG streams
- Enhanced documentation with `@details` section
- Added RNG state restoration

### Documentation Updates
- Explained cross-platform reproducibility approach
- Documented RNG algorithms used (Mersenne-Twister + L'Ecuyer-CMRG)
- Clarified that results are reproducible regardless of thread count

### Backward Compatibility
- Function signature unchanged
- Same behavior for users (just more reproducible)
- RNG state restored after function completes

## Expected Outcomes

1. **Cross-Platform Reproducibility**: Same seed produces identical results on macOS, Linux, and Windows
2. **Thread Independence**: Results are the same regardless of number of threads used
3. **Test Stability**: Vignette workflow tests pass consistently on all platforms
4. **User Confidence**: Scientists can trust that their analyses are reproducible

## References

1. R Documentation: `?RNGkind`, `?clusterSetRNGStream`
2. [ranger issue #533](https://github.com/imbs-hl/ranger/issues/533) - Similar cross-platform RNG issues
3. [Parallel RNG in R](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
4. L'Ecuyer, P. (1999). Good parameters and implementations for combined multiple recursive random number generators. Operations Research, 47, 159-164.
