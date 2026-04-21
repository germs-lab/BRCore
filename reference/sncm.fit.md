# Fit Sloan Neutral Community Model (SNCM)

This function fits the Sloan Neutral Community Model to community data
to assess neutral processes in community assembly. It compares the
neutral model against binomial and Poisson null models using various
statistical measures.

## Usage

``` r
sncm.fit(spp, pool = NULL, stats = TRUE, taxon = NULL)
```

## Arguments

- spp:

  A matrix or data frame where rows represent communities (samples) and
  columns represent species/taxa. Values should be abundances or counts.

- pool:

  Optional. A matrix or data frame representing the source pool for
  calculating relative abundances. If NULL, the source pool is
  calculated from the spp matrix. Default is NULL.

- stats:

  Logical. If TRUE, returns fit statistics including AIC, BIC,
  R-squared, and RMSE for model comparison. If FALSE, returns predicted
  vs observed frequencies with confidence intervals. Default is TRUE.

- taxon:

  Optional. A data frame containing taxonomic information to merge with
  results when stats=FALSE. Should have row names matching species
  names. Default is NULL.

## Value

If stats=TRUE, returns a data frame with model fit statistics including:

- m: Migration rate parameter from NLS fit

- m.ci: Confidence interval for m parameter

- m.mle: Migration rate from maximum likelihood estimation

- maxLL, binoLL, poisLL: Log-likelihood values for SNCM, binomial, and
  Poisson models

- Rsqr, Rsqr.bino, Rsqr.pois: R-squared values for each model

- RMSE, RMSE.bino, RMSE.pois: Root mean squared error for each model

- AIC, BIC: Information criteria for model selection

- N: Average number of individuals per community

- Samples: Number of samples/communities

- Richness: Number of taxa analyzed

- Detect: Detection limit (1/N)

If stats=FALSE, returns a data frame with observed and predicted
frequencies along with confidence intervals for visualization.

## Details

The function implements the Sloan Neutral Community Model which predicts
the occurrence frequency of taxa based on their relative abundance in
the source pool and a migration parameter (m). The model assumes neutral
processes govern community assembly. Three models are compared: SNCM
(beta distribution), binomial null model, and Poisson null model.

## References

Sloan, W.T. et al. (2006) Quantifying the roles of immigration and
chance in shaping prokaryote community structure. Environmental
Microbiology, 8, 732-740.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example community data
spp_data <- matrix(rpois(100 * 50, 5), nrow = 100, ncol = 50)
colnames(spp_data) <- paste0("Species_", 1:50)

# Fit SNCM and get statistics
fit_stats <- sncm.fit(spp_data, stats = TRUE)
print(fit_stats)

# Get predictions for plotting
predictions <- sncm.fit(spp_data, stats = FALSE)
} # }
```
