# Predictive performance

Calculate Bayesian-\\R^2\\, LOO-\\R^2\\, and leave-one-out predictive
metrics for a regression model. Currently this is only calculated for
IPD regression models. Any models with both IPD and AgD (e.g. ML-NMR)
will ignore the AgD part for this calculation.

## Usage

``` r
# S3 method for class 'stan_nma'
bayes_R2(object, ..., probs = c(0.025, 0.5, 0.975), summary = TRUE)

# S3 method for class 'stan_nma'
loo_R2(object, ..., probs = c(0.025, 0.5, 0.975), summary = TRUE)

# S3 method for class 'stan_nma'
loo_predictive_metric(x, metric, ...)
```

## Arguments

- object, x:

  A `stan_nma` object.

- ...:

  Not used.

- probs:

  Numeric vector of quantiles of interest to present in computed
  summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`

- summary:

  Logical, calculate posterior summaries? Default `TRUE`.

- metric:

  The type of predictive metric to be calculated. See
  [`loo::loo_predictive_metric()`](https://mc-stan.org/loo/reference/loo_predictive_metric.html).
  The default is root mean squared error (`"rmse"`) for continuous or
  poisson outcomes, and true positive prediction rate (`"acc"`) for
  binary outcomes.

## Value

For
[`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html)
and
[`loo_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html): a
[nma_summary](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
object if `summary = TRUE`, otherwise a 3D MCMC array of samples of the
Bayesian or LOO \\R^2\\. For `loo_predictive_metric()`: a list with the
estimate and standard error.
