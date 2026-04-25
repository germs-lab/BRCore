# Set up an appropriate parallel backend

Creates a parallel cluster if in an HPC environment or returns the
number of cores to use if on a local machine.

## Usage

``` r
setup_parallel_backend(default = 1)
```

## Arguments

- default:

  Default number of cores to use if detection fails

## Value

Either a cluster object (HPC) or integer number of cores (local)

## Examples

``` r
setup_parallel_backend()
#> ℹ Local environment. Using 3 core(s)/worker(s) via parallelly::availableCores().
#> Keeping 1 core/worker available for system.
#> system 
#>      3 
```
