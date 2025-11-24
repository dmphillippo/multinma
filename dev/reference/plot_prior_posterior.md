# Plot prior vs posterior distribution

Produce plots comparing the prior and posterior distributions of model
parameters.

## Usage

``` r
plot_prior_posterior(
  x,
  ...,
  prior = NULL,
  post_args = list(),
  prior_args = list(),
  overlay = c("prior", "posterior"),
  ref_line = NA_real_
)
```

## Arguments

- x:

  A `stan_nma` object

- ...:

  Additional arguments passed on to methods

- prior:

  Character vector selecting the prior and posterior distribution(s) to
  plot. May include `"intercept"`, `"trt"`, `"het"`, `"reg"`, `"aux"`,
  `"class_mean"` or `"class_sd"` as appropriate.

- post_args:

  List of arguments passed on to
  [ggplot2::geom_histogram](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
  to control plot output for the posterior distribution

- prior_args:

  List of arguments passed on to
  [ggplot2::geom_path](https://ggplot2.tidyverse.org/reference/geom_path.html)
  to control plot output for the prior distribution. Additionally, `n`
  controls the number of points the density curve is evaluated at
  (default `500`), and `p_limits` controls the endpoints of the curve as
  quantiles (default `c(.001, .999)`).

- overlay:

  String, should prior or posterior be shown on top? Default `"prior"`.

- ref_line:

  Numeric vector of positions for reference lines, by default no
  reference lines are drawn

## Value

A `ggplot` object.

## Details

Prior distributions are displayed as lines, posterior distributions are
displayed as histograms.

## Examples

``` r
## Smoking cessation NMA
# \donttest{
# Run smoking RE NMA example if not already available
if (!exists("smk_fit_RE")) example("example_smk_re", run.donttest = TRUE)
# }
# \donttest{
# Plot prior vs. posterior, by default all parameters are plotted
plot_prior_posterior(smk_fit_RE)


# Plot prior vs. posterior for heterogeneity SD only
plot_prior_posterior(smk_fit_RE, prior = "het")


# Customise plot
plot_prior_posterior(smk_fit_RE, prior = "het",
                     prior_args = list(colour = "darkred", linewidth = 2),
                     post_args = list(alpha = 0.6))

# }
```
