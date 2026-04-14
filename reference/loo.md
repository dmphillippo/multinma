# Model comparison using the `loo` package

The [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`waic()`](https://mc-stan.org/loo/reference/waic.html) functions from
the `loo` package may be called directly on
[stan_nma](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
and
[stan_mlnmr](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
objects.

## Usage

``` r
# S3 method for class 'stan_nma'
loo(x, ...)

# S3 method for class 'stan_nma'
waic(x, ...)
```

## Arguments

- x:

  An object of class
  [stan_nma](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
  or
  [stan_mlnmr](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)

- ...:

  Further arguments to `loo()` or
  [`waic()`](https://mc-stan.org/loo/reference/waic.html)
