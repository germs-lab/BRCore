# BRCore Documentation

This directory contains technical documentation for BRCore implementation details and design decisions.

## Contents

### [RNG_REPRODUCIBILITY.md](RNG_REPRODUCIBILITY.md)
Comprehensive documentation on the cross-platform Random Number Generation (RNG) reproducibility solution implemented for `multi_rarefy()` and `identify_core()` functions.

**Topics covered:**
- Problem analysis of cross-platform RNG differences
- Root causes (RNG algorithms, parallel execution, BLAS libraries)
- Implemented solution using explicit RNG algorithms
- Alternative approaches considered
- Testing strategy
- References and related issues

**Key audience:** 
- Package maintainers
- Contributors working on random sampling functions
- Users experiencing reproducibility issues
- Researchers needing to understand the technical approach

## Contributing

When adding new documentation:
1. Create a descriptive filename using SCREAMING_SNAKE_CASE
2. Include a clear problem statement and solution
3. Document alternatives considered
4. Provide references and examples
5. Update this README with a link to the new document
