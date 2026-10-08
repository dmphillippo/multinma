# Posterior summaries from `stan_nma` objects

Posterior summaries of model parameters in `stan_nma` objects may be
produced using the [`summary()`](https://rdrr.io/r/base/summary.html)
method and plotted with the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method. NOTE:
To produce relative effects, absolute predictions, or posterior ranks,
see
[`relative_effects()`](https://dmphillippo.github.io/multinma/dev/reference/relative_effects.md),
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md),
[`posterior_ranks()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md),
[`posterior_rank_probs()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md).

## Usage

``` r
# S3 method for class 'stan_nma'
summary(object, ..., pars, include, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# S3 method for class 'stan_nma'
plot(
  x,
  ...,
  pars,
  include,
  stat = "pointinterval",
  orientation = c("horizontal", "vertical", "y", "x"),
  ref_line = NA_real_
)
```

## Arguments

- ...:

  Additional arguments passed on to other methods

- pars, include:

  See
  [`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)

- probs:

  Numeric vector of specifying quantiles of interest, default
  `c(0.025, 0.25, 0.5, 0.75, 0.975)`

- x, object:

  A `stan_nma` object

- stat:

  Character string specifying the `ggdist` plot stat to use, default
  `"pointinterval"`

- orientation:

  Whether the `ggdist` geom is drawn horizontally (`"horizontal"`) or
  vertically (`"vertical"`), default `"horizontal"`

- ref_line:

  Numeric vector of positions for reference lines, by default no
  reference lines are drawn

- summary:

  Logical, calculate posterior summaries? Default `TRUE`.

## Value

A
[nma_summary](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
object

## Details

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method is a
shortcut for `plot(summary(stan_nma))`. For details of plotting options,
see
[`plot.nma_summary()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md).

## See also

[`plot.nma_summary()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md),
[`relative_effects()`](https://dmphillippo.github.io/multinma/dev/reference/relative_effects.md),
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md),
[`posterior_ranks()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md),
[`posterior_rank_probs()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md)

## Examples

``` r
## Smoking cessation
# \donttest{
# Run smoking RE NMA example if not already available
if (!exists("smk_fit_RE")) example("example_smk_re", run.donttest = TRUE)
# }
# \donttest{
# Summary and plot of all model parameters
summary(smk_fit_RE)
#>                                    mean   sd  2.5%   25%   50%   75% 97.5%
#> mu[1]                             -2.78 0.33 -3.44 -2.99 -2.77 -2.56 -2.17
#> mu[2]                             -2.54 0.78 -4.12 -3.02 -2.54 -2.03 -1.04
#> mu[3]                             -2.14 0.12 -2.39 -2.22 -2.14 -2.06 -1.90
#> mu[4]                             -4.05 0.56 -5.20 -4.41 -4.03 -3.66 -3.03
#> mu[5]                             -2.15 0.14 -2.43 -2.25 -2.15 -2.05 -1.89
#> mu[6]                             -3.41 0.71 -4.92 -3.84 -3.35 -2.91 -2.18
#> mu[7]                             -3.02 0.45 -3.98 -3.30 -3.00 -2.71 -2.23
#> mu[8]                             -2.70 0.61 -4.04 -3.06 -2.65 -2.27 -1.64
#> mu[9]                             -1.84 0.42 -2.70 -2.11 -1.82 -1.55 -1.07
#> mu[10]                            -2.08 0.12 -2.32 -2.16 -2.08 -2.00 -1.85
#> mu[11]                            -3.63 0.23 -4.11 -3.78 -3.62 -3.46 -3.19
#> mu[12]                            -2.22 0.13 -2.48 -2.31 -2.22 -2.13 -1.98
#> mu[13]                            -2.68 0.45 -3.61 -2.97 -2.66 -2.37 -1.84
#> mu[14]                            -2.41 0.23 -2.89 -2.56 -2.40 -2.25 -1.98
#> mu[15]                            -2.69 0.75 -4.30 -3.15 -2.63 -2.17 -1.39
#> mu[16]                            -2.62 0.34 -3.32 -2.83 -2.60 -2.38 -2.00
#> mu[17]                            -2.38 0.11 -2.60 -2.45 -2.37 -2.30 -2.17
#> mu[18]                            -2.57 0.27 -3.10 -2.74 -2.56 -2.39 -2.07
#> mu[19]                            -1.90 0.12 -2.14 -1.98 -1.90 -1.82 -1.66
#> mu[20]                            -2.80 0.13 -3.05 -2.88 -2.80 -2.72 -2.56
#> mu[21]                            -1.13 0.81 -2.77 -1.64 -1.14 -0.59  0.45
#> mu[22]                            -2.41 0.85 -4.13 -2.96 -2.38 -1.86 -0.78
#> mu[23]                            -2.31 0.83 -3.98 -2.84 -2.31 -1.78 -0.68
#> mu[24]                            -2.80 0.86 -4.56 -3.35 -2.81 -2.23 -1.12
#> d[Group counselling]               1.10 0.44  0.27  0.81  1.09  1.38  2.03
#> d[Individual counselling]          0.84 0.25  0.37  0.67  0.84  1.00  1.35
#> d[Self-help]                       0.50 0.41 -0.30  0.24  0.49  0.76  1.30
#> tau                                0.84 0.19  0.54  0.71  0.81  0.94  1.27
#> delta[1: Individual counselling]   1.07 0.38  0.33  0.81  1.07  1.33  1.83
#> delta[1: Group counselling]        0.37 0.42 -0.44  0.09  0.37  0.65  1.19
#> delta[2: Self-help]                0.64 0.80 -0.92  0.13  0.62  1.13  2.27
#> delta[2: Individual counselling]   0.73 0.79 -0.80  0.22  0.72  1.23  2.37
#> delta[2: Group counselling]        0.96 0.79 -0.56  0.46  0.97  1.47  2.52
#> delta[3: Individual counselling]   2.16 0.14  1.88  2.07  2.16  2.26  2.45
#> delta[4: Individual counselling]   0.90 0.58 -0.20  0.51  0.89  1.29  2.06
#> delta[5: Individual counselling]   0.44 0.16  0.14  0.33  0.43  0.54  0.75
#> delta[6: Individual counselling]   1.72 0.72  0.43  1.22  1.67  2.17  3.22
#> delta[7: Individual counselling]   2.15 0.48  1.27  1.80  2.13  2.45  3.18
#> delta[8: Individual counselling]   1.64 0.62  0.55  1.22  1.60  2.02  3.02
#> delta[9: Individual counselling]   0.59 0.46 -0.30  0.28  0.59  0.89  1.53
#> delta[10: Self-help]               0.00 0.17 -0.32 -0.11  0.01  0.12  0.32
#> delta[11: Self-help]               0.41 0.30 -0.17  0.22  0.41  0.61  1.00
#> delta[12: Individual counselling]  0.41 0.17  0.09  0.30  0.42  0.53  0.75
#> delta[13: Individual counselling]  0.40 0.51 -0.57  0.06  0.40  0.74  1.41
#> delta[14: Individual counselling]  0.63 0.29  0.07  0.44  0.63  0.82  1.20
#> delta[15: Group counselling]       2.15 0.78  0.76  1.61  2.10  2.61  3.84
#> delta[16: Self-help]               0.66 0.40 -0.10  0.39  0.65  0.92  1.46
#> delta[17: Individual counselling]  0.55 0.14  0.29  0.46  0.55  0.64  0.83
#> delta[18: Individual counselling]  0.03 0.31 -0.57 -0.18  0.02  0.23  0.65
#> delta[19: Individual counselling] -0.19 0.17 -0.53 -0.31 -0.19 -0.08  0.14
#> delta[20: Individual counselling]  0.08 0.19 -0.29 -0.04  0.08  0.21  0.45
#> delta[21: Self-help]               0.70 0.81 -0.90  0.18  0.69  1.21  2.30
#> delta[21: Individual counselling]  0.65 0.80 -0.93  0.12  0.66  1.16  2.24
#> delta[22: Self-help]               0.32 0.84 -1.31 -0.24  0.31  0.85  2.02
#> delta[22: Group counselling]       1.29 0.85 -0.32  0.73  1.27  1.83  2.99
#> delta[23: Individual counselling]  0.65 0.81 -0.92  0.13  0.64  1.17  2.24
#> delta[23: Group counselling]       1.26 0.84 -0.33  0.72  1.26  1.78  2.95
#> delta[24: Individual counselling]  1.04 0.83 -0.59  0.50  1.04  1.56  2.65
#> delta[24: Group counselling]       0.90 0.87 -0.85  0.34  0.91  1.46  2.56
#>                                   Bulk_ESS Tail_ESS Rhat
#> mu[1]                                 5035     2967 1.00
#> mu[2]                                 2639     2638 1.00
#> mu[3]                                10078     2794 1.00
#> mu[4]                                 4478     2698 1.00
#> mu[5]                                 7861     2879 1.00
#> mu[6]                                 3799     2249 1.00
#> mu[7]                                 4114     2659 1.00
#> mu[8]                                 3843     2423 1.00
#> mu[9]                                 5106     2848 1.00
#> mu[10]                                8825     2537 1.00
#> mu[11]                                7680     3062 1.00
#> mu[12]                                7624     2502 1.00
#> mu[13]                                5320     3097 1.00
#> mu[14]                                6093     2503 1.00
#> mu[15]                                3369     2196 1.00
#> mu[16]                                5719     2751 1.00
#> mu[17]                                7684     2965 1.00
#> mu[18]                                6324     3093 1.00
#> mu[19]                                8218     3036 1.00
#> mu[20]                                9148     3163 1.00
#> mu[21]                                2842     2078 1.00
#> mu[22]                                2937     2308 1.00
#> mu[23]                                2967     2934 1.00
#> mu[24]                                3002     2605 1.00
#> d[Group counselling]                  1924     2375 1.00
#> d[Individual counselling]             1147     1549 1.00
#> d[Self-help]                          1830     2353 1.00
#> tau                                   1154     1750 1.01
#> delta[1: Individual counselling]      5317     3164 1.00
#> delta[1: Group counselling]           4364     3214 1.00
#> delta[2: Self-help]                   2670     2439 1.00
#> delta[2: Individual counselling]      2763     2699 1.00
#> delta[2: Group counselling]           2675     2718 1.00
#> delta[3: Individual counselling]      7520     3320 1.00
#> delta[4: Individual counselling]      4468     3512 1.00
#> delta[5: Individual counselling]      6777     3243 1.00
#> delta[6: Individual counselling]      3705     2598 1.00
#> delta[7: Individual counselling]      3847     2991 1.00
#> delta[8: Individual counselling]      3911     2620 1.00
#> delta[9: Individual counselling]      4873     3091 1.00
#> delta[10: Self-help]                  5773     3433 1.00
#> delta[11: Self-help]                  5952     3386 1.00
#> delta[12: Individual counselling]     5655     2758 1.00
#> delta[13: Individual counselling]     4742     3211 1.00
#> delta[14: Individual counselling]     5664     2574 1.00
#> delta[15: Group counselling]          3197     2709 1.00
#> delta[16: Self-help]                  5565     2871 1.00
#> delta[17: Individual counselling]     6941     3331 1.00
#> delta[18: Individual counselling]     5975     3210 1.00
#> delta[19: Individual counselling]     5411     3525 1.00
#> delta[20: Individual counselling]     6023     3523 1.00
#> delta[21: Self-help]                  2822     2125 1.00
#> delta[21: Individual counselling]     2852     2257 1.00
#> delta[22: Self-help]                  3046     2574 1.00
#> delta[22: Group counselling]          2882     2083 1.00
#> delta[23: Individual counselling]     2998     2787 1.00
#> delta[23: Group counselling]          3023     2913 1.00
#> delta[24: Individual counselling]     3072     2753 1.00
#> delta[24: Group counselling]          3106     2752 1.00
plot(smk_fit_RE)


# Summary and plot of heterogeneity tau only
summary(smk_fit_RE, pars = "tau")
#>     mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> tau 0.84 0.19 0.54 0.71 0.81 0.94  1.27     1154     1750 1.01
plot(smk_fit_RE, pars = "tau")


# Customising plot output
plot(smk_fit_RE,
     pars = c("d", "tau"),
     stat = "halfeye",
     ref_line = 0)

# }
```
