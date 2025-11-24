# Working with 3D MCMC arrays

3D MCMC arrays (Iterations, Chains, Parameters) are produced by
[`as.array()`](https://rdrr.io/r/base/array.html) methods applied to
`stan_nma` or `nma_summary` objects.

## Usage

``` r
# S3 method for class 'mcmc_array'
summary(object, ..., probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# S3 method for class 'mcmc_array'
print(x, ...)

# S3 method for class 'mcmc_array'
plot(x, ...)

# S3 method for class 'mcmc_array'
names(x)

# S3 method for class 'mcmc_array'
names(x) <- value
```

## Arguments

- ...:

  Further arguments passed to other methods

- probs:

  Numeric vector of quantiles of interest

- x, object:

  A 3D MCMC array of class `mcmc_array`

- value:

  Character vector of replacement parameter names

## Value

The [`summary()`](https://rdrr.io/r/base/summary.html) method returns a
[nma_summary](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
object, the [`print()`](https://rdrr.io/r/base/print.html) method
returns `x` invisibly. The
[`names()`](https://rdrr.io/r/base/names.html) method returns a
character vector of parameter names, and `names()<-` returns the object
with updated parameter names. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method is a
shortcut for `plot(summary(x), ...)`, passing all arguments on to
[`plot.nma_summary()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md).

## Examples

``` r
## Smoking cessation
# \donttest{
# Run smoking RE NMA example if not already available
if (!exists("smk_fit_RE")) example("example_smk_re", run.donttest = TRUE)
# }
# \donttest{
# Working with arrays of posterior draws (as mcmc_array objects) is
# convenient when transforming parameters

# Transforming log odds ratios to odds ratios
LOR_array <- as.array(relative_effects(smk_fit_RE))
OR_array <- exp(LOR_array)

# mcmc_array objects can be summarised to produce a nma_summary object
smk_OR_RE <- summary(OR_array)

# This can then be printed or plotted
smk_OR_RE
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> d[Group counselling]      3.36 1.71 1.29 2.28 3.01 4.02  7.79     1541     2173
#> d[Individual counselling] 2.42 0.62 1.48 1.98 2.32 2.74  3.91      985     1773
#> d[Self-help]              1.77 0.80 0.74 1.25 1.62 2.10  3.74     1483     1676
#>                           Rhat
#> d[Group counselling]         1
#> d[Individual counselling]    1
#> d[Self-help]                 1
plot(smk_OR_RE, ref_line = 1)


# Transforming heterogeneity SD to variance
tau_array <- as.array(smk_fit_RE, pars = "tau")
tausq_array <- tau_array^2

# Correct parameter names
names(tausq_array) <- "tausq"

# Summarise
summary(tausq_array)
#>       mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> tausq 0.75 0.35  0.3 0.51 0.67 0.91  1.62     1116     1643    1
# }
```
