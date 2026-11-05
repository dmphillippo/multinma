# Baseline synthesis models

Runs a baseline synthesis model, placing a random effect on the
study-specific intercepts to obtain a pooled estimate of the absolute
outcomes on the network reference treatment (or the baseline treatment
specified by `baseline_trt`). This can be combined with a regression
model to obtain an adjusted estimate. Estimates of both the overall mean
and standard deviation (`baseline_mean` and `baseline_sd`) as well as
the predictive distribution (`baseline_new`) are returned.

## Usage

``` r
baseline_synthesis(
  network,
  baseline_trt = NULL,
  prior_intercept = .default(normal(scale = 100)),
  prior_intercept_sd = .default(half_normal(scale = 5)),
  ...
)
```

## Arguments

- network:

  A `multinma` network object.

- baseline_trt:

  Baseline treatment on which to place the random effect for baseline
  synthesis. If not specified, the network reference treatment will be
  used.

- prior_intercept:

  Prior for the mean of the random effects distribution on the
  study-specific intercepts.

- prior_intercept_sd:

  Prior for standard deviation of the random effects distribution on the
  study-specific intercepts.

- ...:

  Additional arguments passed directly to
  [`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md).

## Value

A `stan_baseline` object, inheriting from
[stan_nma](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md).

## Details

The baseline synthesis model places a random effect on the
study-specific intercept parameters of the NMA or ML-NMR model fitted by
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md):
\$\$\mu_j \sim \mathrm{N}(m, \tau\_\mu)\$\$ The prior distribution on
the mean \\m\\ is specified by the `prior_intercept` argument, and the
prior for \\\tau\_\mu\\ by `prior_intercept_sd`.
