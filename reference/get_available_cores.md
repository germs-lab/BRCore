# Get the number of available CPU cores based on environment

Determines the appropriate number of cores to use based on whether we're
in an HPC environment or a local machine.

## Usage

``` r
get_available_cores(default = 1)
```

## Arguments

- default:

  Default number of cores to use if detection fails

## Value

Integer number of cores available

## Examples

``` r
get_available_cores()
#> ℹ Local environment. Using 3 core(s)/worker(s) via parallelly::availableCores().
#> Keeping 1 core/worker available for system.
#> system 
#>      3 
get_available_cores(default = 2)
#> ℹ Local environment. Using 3 core(s)/worker(s) via parallelly::availableCores().
#> Keeping 1 core/worker available for system.
#> system 
#>      3 
```
