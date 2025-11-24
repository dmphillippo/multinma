# Convert samples into arrays, matrices, or data frames

Samples (post warm-up) from a `stan_nma` model object can be coerced
into an array, matrix, or data frame.

## Usage

``` r
# S3 method for class 'stan_nma'
as.array(x, ..., pars, include = TRUE)

# S3 method for class 'stan_nma'
as.data.frame(x, ..., pars, include = TRUE)

# S3 method for class 'stan_nma'
as_tibble(x, ..., pars, include = TRUE)

# S3 method for class 'stan_nma'
as.tibble(x, ..., pars, include = TRUE)

# S3 method for class 'stan_nma'
as.matrix(x, ..., pars, include = TRUE)
```

## Arguments

- x:

  A `stan_nma` object

- ...:

  Additional arguments passed to
  [`as.array.stanfit()`](https://mc-stan.org/rstan/reference/stanfit2array-method.html)

- pars:

  Optional character vector of parameter names to include in output. If
  not specified, all parameters are used.

- include:

  Logical, are parameters in `pars` to be included (`TRUE`, default) or
  excluded (`FALSE`)?

## Value

The [`as.array()`](https://rdrr.io/r/base/array.html) method produces a
3D array \[Iteration, Chain, Parameter\] containing posterior samples of
each parameter (as class
[mcmc_array](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)).
This has the side effect of enabling
[`bayesplot`](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
functions to seamlessly work on `stan_nma` objects.

The [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
method produces a data frame containing posterior samples of each
parameter, combined over all chains.

The [`as.matrix()`](https://rdrr.io/r/base/matrix.html) method produces
a matrix containing posterior samples of each parameter, combined over
all chains.
