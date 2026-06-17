# Bind chains

Methods to combine chains from fitted
[stan_nma](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
model objects, or
[mcmc_array](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
3D MCMC arrays. `cbind.stan_nma()` is an alias for `bind_chains()`.

## Usage

``` r
bind_chains(...)

# S3 method for class 'stan_nma'
cbind(...)

# S3 method for class 'mcmc_array'
cbind(...)
```

## Arguments

- ...:

  Multiple fitted
  [stan_nma](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
  model objects (for `bind_chains()`) or
  [mcmc_array](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  arrays (for `cbind.mcmc_array()`). `bind_chains()` also accepts a
  single list as input, containing the objects to combine.

## Value

A
[stan_nma](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
or
[mcmc_array](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
object.

## Details

Objects to combine must have the same number of iterations, warmup, and
thinning, and all parameter names must match. For `bind_chains()`, this
means that both the model specification and input data must be the same.
