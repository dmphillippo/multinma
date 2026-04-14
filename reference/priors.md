# Prior distributions

These functions are used to specify prior distributions for the model
parameters.

## Usage

``` r
normal(location = 0, scale)

half_normal(scale)

log_normal(location, scale)

cauchy(location = 0, scale)

half_cauchy(scale)

student_t(location = 0, scale, df)

half_student_t(scale, df)

log_student_t(location, scale, df)

exponential(scale = 1/rate, rate = 1/scale)

flat()
```

## Arguments

- location:

  Prior location. Typically prior mean (see details).

- scale:

  Prior scale. Typically prior standard deviation (see details).

- df:

  Prior degrees of freedom.

- rate:

  Prior rate.

## Value

Object of class
[nma_prior](https://dmphillippo.github.io/multinma/reference/nma_prior-class.md).

## Details

The `location` and `scale` parameters are typically the prior mean and
standard deviation, with the following exceptions:

- For the Cauchy distribution `location` is the prior median and `scale`
  is the prior scale.

- For the log-Normal distribution, `location` and `scale` are the prior
  mean and standard deviation of the logarithm.

### Compatibility with model parameters

The following table summarises which prior distributions may be used
with which model parameters. Essentially, priors that take only
non-negative values (e.g. half-Normal) may only be used for non-negative
parameters (heterogeneity SD/variance/precision, and any auxiliary
parameter). If a real-valued prior distribution is specified for a
non-negative parameter, it will be truncated at 0 to be non-negative.

|                                       |                                 |                                   |                               |                                         |                                     |
|---------------------------------------|---------------------------------|-----------------------------------|-------------------------------|-----------------------------------------|-------------------------------------|
|                                       | **Intercept** `prior_intercept` | **Treatment effects** `prior_trt` | **Heterogeneity** `prior_het` | **Regression coefficients** `prior_reg` | **Auxiliary parameter** `prior_aux` |
| **Normal** `normal()`                 | Yes                             | Yes                               | Yes                           | Yes                                     | Yes                                 |
| **half-Normal** `half_normal()`       | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **log-Normal** `log_normal()`         | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **Cauchy** `cauchy()`                 | Yes                             | Yes                               | Yes                           | Yes                                     | Yes                                 |
| **half-Cauchy** `half_cauchy()`       | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **Student t** `student_t()`           | Yes                             | Yes                               | Yes                           | Yes                                     | Yes                                 |
| **half-Student t** `half_student_t()` | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **log-Student t** `log_student_t()`   | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **Exponential** `exponential()`       | \-                              | \-                                | Yes                           | \-                                      | Yes                                 |
| **Flat** `flat()`                     | Yes                             | Yes                               | Yes                           | Yes                                     | Yes                                 |

The `flat()` prior is a special case where no prior information is added
to the model, resulting in an implicit flat uniform prior distribution
over the entire support for a parameter. This will be an improper prior
if the parameter is unbounded, and is not generally advised. See the
[Stan user's
guide](https://mc-stan.org/docs/stan-users-guide/some-differences-in-the-statistical-models-that-are-allowed.html)
for more details.

## See also

[`summary.nma_prior()`](https://dmphillippo.github.io/multinma/reference/summary.nma_prior.md)
for summarising details of prior distributions.
[`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/reference/plot_prior_posterior.md)
for plots comparing the prior and posterior distributions of model
parameters.
