# Matrix of plots for a `stan_nma` object

A [`pairs()`](https://rdrr.io/r/graphics/pairs.html) method for
`stan_nma` objects, which calls
[`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html)
on the underlying `stanfit` object.

## Usage

``` r
# S3 method for class 'stan_nma'
pairs(x, ..., pars, include = TRUE)
```

## Arguments

- x:

  An object of class `stan_nma`

- ...:

  Other arguments passed to
  [`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html)

- pars:

  Optional character vector of parameter names to include in output. If
  not specified, all parameters are used.

- include:

  Logical, are parameters in `pars` to be included (`TRUE`, default) or
  excluded (`FALSE`)?

## Value

A grid of ggplot objects produced by
[`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html).

## Examples

``` r
if (FALSE) { # \dontrun{
## Parkinson's mean off time reduction
park_net <- set_agd_arm(parkinsons,
                        study = studyn,
                        trt = trtn,
                        y = y,
                        se = se,
                        sample_size = n)

# Fitting a RE model
park_fit_RE <- nma(park_net,
                   trt_effects = "random",
                   prior_intercept = normal(scale = 100),
                   prior_trt = normal(scale = 100),
                   prior_het = half_normal(scale = 5))

# We see a small number of divergent transition errors
# These do not go away entirely when adapt_delta is increased

# Try to diagnose with a pairs plot
pairs(park_fit_RE, pars = c("mu[4]", "d[3]", "delta[4: 3]", "tau"))

# Transforming tau onto log scale
pairs(park_fit_RE, pars = c("mu[4]", "d[3]", "delta[4: 3]", "tau"),
      transformations = list(tau = "log"))

# The divergent transitions occur in the upper tail of the heterogeneity
# standard deviation. In this case, with only a small number of studies, there
# is not very much information to estimate the heterogeneity standard deviation
# and the prior distribution may be too heavy-tailed. We could consider a more
# informative prior distribution for the heterogeneity variance to aid
# estimation.
} # }
```
