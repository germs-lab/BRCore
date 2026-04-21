# Detect HPC environment type

Detects if code is running in an HPC environment and identifies which
scheduler.

## Usage

``` r
detect_hpc_type()
```

## Value

Character string identifying the HPC environment type: "slurm", "pbs",
"sge", "lsf", "condor", "generic_hpc", or "local"

## Examples

``` r
detect_hpc_type()
#> [1] "local"
```
