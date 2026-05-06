# Fit Sloan Neutral Community Model (SNCM)

The function implements the Sloan Neutral Community Model which predicts
the occurrence frequency of taxa based on their relative abundance in
the source pool and a migration parameter (`m`). It compares the neutral
model against binomial and Poisson null models using various statistical
measures.

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
  calculated from the spp matrix. Default is `NULL`.

- stats:

  Logical. If `TRUE`, returns fit statistics including AIC, BIC,
  R-squared, and RMSE for model comparison. If `FALSE`, returns
  predicted vs. observed frequencies with confidence intervals. Default
  is `TRUE`.

- taxon:

  Optional. A data frame containing taxonomic information to merge with
  results when `stats = FALSE`. Should have row names matching species
  names. Default is `NULL`.

## Value

If `stats = TRUE`, returns a data frame with model fit statistics
including:

- `m`: Migration rate parameter from NLS fit

- `m.ci`: Confidence interval for m parameter

- `m.mle`: Migration rate from maximum likelihood estimation

- `maxLL`, `binoLL`, `poisLL`: Log-likelihood values for SNCM, binomial,
  and Poisson models

- `Rsqr`, `Rsqr.bino`, `Rsqr.pois`: R-squared values for each model

- `RMSE`, `RMSE.bino`, `RMSE.pois`: Root mean squared error for each
  model

- `AIC`, `BIC`: Information criteria for model selection

- `N`: Average number of individuals per community

- `Samples`: Number of samples/communities

- `Richness`: Number of taxa analyzed

- `Detect`: Detection limit (1/N)

If `stats = FALSE`, returns a data frame with observed and predicted
frequencies along with confidence intervals for visualization.

## Details

The model assumes neutral processes govern community assembly. Three
models are compared: SNCM (beta distribution), binomial null model, and
Poisson null model.

## References

Sloan WT, Lunn M, Woodcock S, Head IM, Nee S, Curtis TP. (2006)
Quantifying the roles of immigration and chance in shaping prokaryote
community structure. Environ Microbiol. 8(4):732-40.
<doi:10.1111/j.1462-2920.2005.00956.x>

## See also

[`plot_neutral_model()`](https://www.germslab.org/BRCore/reference/plot_neutral_model.md),
[`fit_neutral_model()`](https://www.germslab.org/BRCore/reference/fit_neutral_model.md)

## Examples

``` r
# Generate community data
set.seed(42)
n_samples <- 100
n_species <- 80

# Log-normal abundances
mean_abund <- rlnorm(n_species, meanlog = 2, sdlog = 1.5)
# Simulate community matrix with multinomial sampling
spp_data <- t(sapply(seq_len(n_samples), function(i) {
  rmultinom(
    1,
    size = sample(500:2000, 1),
    prob = mean_abund / sum(mean_abund)
  )[, 1]
}))
colnames(spp_data) <- paste0("Species_", seq_len(n_species))

fit_stats <- sncm.fit(spp_data, stats = TRUE)
#> Waiting for profiling to be done...
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> ℹ Neutral model fitting:
#> • Average individuals per community (N): 1264.75
#> • Binomial model using rounded N: 1265
#> • Poisson model using N: 1264.75
#> • Maximum likelihood estimation using N: 1264.75, and starting parameters: mu =
#>   0, sigma = 0.1
predictions <- sncm.fit(spp_data, stats = FALSE)
#> Waiting for profiling to be done...
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> ℹ Neutral model fitting:
#> • Average individuals per community (N): 1264.75
#> • Binomial model using rounded N: 1265
#> • Poisson model using N: 1264.75
#> • Maximum likelihood estimation using N: 1264.75, and starting parameters: mu =
#>   0, sigma = 0.1

```
