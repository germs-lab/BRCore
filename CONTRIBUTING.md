# Contributing to BRCore

Thank you for your interest in contributing to BRCore! This document
provides general guidelines for contributing to the package.

## Table of Contents

- [Getting Help](#getting-help)
- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Workflow](#development-workflow)
- [Additional Resources](#additional-resources)

## Getting Help

- **Questions:** Open a discussion on GitHub Discussions

- **Bug reports:** Open an issue on [GitHub
  Issues](https://github.com/germs-lab/BRCore/issues)

- **Contact:** Reach out to package maintainers listed in `DESCRIPTION`

## Code of Conduct

BRCore is committed to providing a welcoming and harassment-free
experience for everyone. We expect all contributors to adhere to our
Code of Conduct is adapted from the Contributor Covenant. Please be
respectful and constructive in all interactions.

## Getting Started

### Prerequisites

- R ≥ 4..0
- Git
- RStudio or Positron (recommended)
- Required packages: `devtools`, `testthat`, `roxygen2`, `pkgdown`

### Setting Up Your Development Environment

1.  **Fork and clone the repository:**

    CLassic way:

    ``` bash
    git clone https://github.com/germs-lab/BRCore.git
    cd BRCore
    ```

You can also use `usethis` R package. See [Pull Request
process](https://tidyverse.tidyverse.org/CONTRIBUTING.html#pull-request-process)
for more information.

2\. **Install development dependencies:**

Install required packages

`{r} install.packages(c("devtools", "testthat", "roxygen2", "pkgdown", "usethis"))`

Install package dependencies

`{r} devtools::install_deps(dependencies = TRUE)`

Load the package for development:

`{r} devtools::load_all()`

## How to Contribute

### Types of Contributions

We welcome several types of contributions:

- **Bug reports:** Found a bug? Let us know! Submit a pull request to
  fix existing issues
- **Feature requests:** Have an idea for a new feature? Implement
  requested features or propose new ones
- **Documentation**: Improve examples, vignettes, or function
  documentation
- **Tests:** Add or improve test coverage

### Reporting Bugs

Before creating a bug report, please:

1.  Check existing issues to avoid duplicates
2.  Use the latest development version to ensure the bug hasn’t been
    fixed
3.  Create a minimal reproducible example (reprex)

When filing a bug report, include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected vs. actual behavior
- Your R session info
  ([`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html))
- A minimal reproducible example using the reprex package

Example: library(reprex)

\`\`\`{r} reprex({ library(BRCore) data(“bcse”)

\# Code that demonstrates the bug multi_rarefy(bcse, depth_level = 1000)

sessionInfo() })


    ### Suggesting Features

    Feature requests are welcome! Please:

    -   **Check existing issues** to see if it's already proposed
    -   **Describe the use case:**, why is this feature needed?
    -   **Provide examples** of how the feature would be used
    -   **Consider backward compatibility:** will this break existing code?

    ## Development Workflow {#development-workflow}

    We follow a Git Flow workflow with `main` and `develop` branches.

    ### Branch Strategy

    -   `main`: Stable releases only

    -   `develop`: Integration branch for features

    -   `feature/*`: New features

    -   `patch/*`: Bug/hot fixes

    -   `release/*`: Release preparation

    ### Creating a Feature Branch

    ``` r
    # Start from develop
    git checkout develop
    git pull origin develop

    # Create a feature branch
    git checkout -b feature/your-feature-name

### Making Changes

1.  **Write code following the tidyverse style guide** (see Additional
    Resources)

2.  **Add tests** for new functionality

3.  **Update documentation** if you change function behavior

4.  **Run checks** before committing

### Committing Changes

As of v1.0.2, we follow a mix of [Semantic
Versioning](https://semver.org/) principles and [Conventional
Commits](https://www.conventionalcommits.org/en/v1.0.0/) format:

    <type>(<scope>): <subject>

    [optional body]

    [optional footer]

#### Types:

- `feat`: New feature (minor version bump)

- `fix`: Bug fix (patch version bump)

- `docs`: Documentation changes test: Adding or updating tests

- `refactor`: Code refactoring without behavior changes

- `style`: Code style changes (formatting)

- `chore`: Maintenance tasks (dependencies, build, CI)

- `perf`: Performance improvements

- `BREAKING CHANGE`: Breaking changes (major version bump)

Examples:

``` r
# New feature
git commit -m "feat(multi_rarefy): add reproducible iteration seeds"

# Bug fix
git commit -m "fix(identify_core): add floating-point tolerance in validation"

# Documentation
git commit -m "docs(vignette): update rarefaction workflow examples"

# Breaking change
git commit -m "feat(identify_core): change return type to list

BREAKING CHANGE: identify_core() now returns a list instead of a data frame"
```

### Submitting Changes

1.  **Ensure your code passes all checks**:

    ``` r

    devtools::document()
    devtools::test()
    devtools::check()
    ```

2.  **Push your branch and create a Pull Request**:

- Base: `develop` (not `main`)

- Compare: `feature/your-feature-name`

- Use a clear, descriptive title

- Reference related issues (e.g., “Closes \#123”)

## Additional Resources

For detailed guidance on coding standards, testing, documentation, and
code review:

- **Contributing to tidyverse packages:**
  <https://tidyverse.tidyverse.org/CONTRIBUTING.html>

- **Tidyverse style guide:** <https://style.tidyverse.org/>

- **Code review principles:** <https://code-review.tidyverse.org/>

## Recognition

All contributors will be acknowledged in:

- Package DESCRIPTION file (`ctb` role)

- Release notes

- Package website (contributors page)
