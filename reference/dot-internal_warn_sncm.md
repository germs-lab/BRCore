# Internal warning handler for sncm.fit likelihood functions

Wraps a function call with a
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) that
captures warnings and re-emits them via
[`cli::cli_warn()`](https://cli.r-lib.org/reference/cli_abort.html).
Deduplicates repeated warnings across optimizer iterations using an
external flag and setter.

## Usage

``` r
.internal_warn_sncm(fun, message = NULL, warned = FALSE, warn_setter = NULL)
```

## Arguments

- fun:

  A zero-argument function (lambda) to evaluate.

- message:

  A string prefix for the warning message.

- warned:

  Logical flag indicating whether a warning has already fired.

- warn_setter:

  A zero-argument function that sets the external `warned` flag to
  `TRUE` via `<<-`.

## Value

The result of `fun()`.
