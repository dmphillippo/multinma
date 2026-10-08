# Marginal treatment effects

Generate population-average marginal treatment effects. These are formed
from population-average absolute predictions, so this function is a
wrapper around
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md).

## Usage

``` r
marginal_effects(
  object,
  ...,
  mtype = c("difference", "ratio", "link"),
  all_contrasts = FALSE,
  trt_ref = NULL,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  predictive_distribution = FALSE,
  summary = TRUE
)
```

## Arguments

- object:

  A `stan_nma` object created by
  [`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md).

- ...:

  Arguments passed to
  [`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md),
  for example to specify the covariate distribution and baseline risk
  for a target population, e.g. `newdata`, `baseline`, and related
  arguments. For survival outcomes, `type` can also be specified to
  determine the quantity from which to form a marginal effect. For
  example, `type = "hazard"` with `mtype = "ratio"` produces marginal
  hazard ratios, `type = "median"` with `mtype = "difference"` produces
  marginal median survival time differences, and so on.

- mtype:

  The type of marginal effect to construct from the average absolute
  effects, either `"difference"` (the default) for a difference of
  absolute effects such as a risk difference, `"ratio"` for a ratio of
  absolute effects such as a risk ratio, or `"link"` for a difference on
  the scale of the link function used in fitting the model such as a
  marginal log odds ratio.

- all_contrasts:

  Logical, generate estimates for all contrasts (`TRUE`), or just the
  "basic" contrasts against the network reference treatment (`FALSE`)?
  Default `FALSE`.

- trt_ref:

  Reference treatment to construct relative effects against, if
  `all_contrasts = FALSE`. By default, relative effects will be against
  the network reference treatment. Coerced to character string.

- probs:

  Numeric vector of quantiles of interest to present in computed
  summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`

- predictive_distribution:

  Logical, when a random effects model has been fitted, should the
  predictive distribution for marginal effects in a new study be
  returned? Default `FALSE`.

- summary:

  Logical, calculate posterior summaries? Default `TRUE`.

## Value

A
[nma_summary](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
object if `summary = TRUE`, otherwise a list containing a 3D MCMC array
of samples and (for regression models) a data frame of study
information.

## Examples

``` r
## Smoking cessation
# \donttest{
# Run smoking RE NMA example if not already available
if (!exists("smk_fit_RE")) example("example_smk_re", run.donttest = TRUE)
# }
# \donttest{
# Marginal risk difference in each study population in the network
marginal_effects(smk_fit_RE, mtype = "difference")
#> ---------------------------------------------------------------------- Study: 1 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[1: Group counselling]      0.11 0.07  0.02 0.06 0.10 0.14  0.27     2129
#> marg[1: Individual counselling] 0.07 0.03  0.02 0.05 0.07 0.09  0.15     1503
#> marg[1: Self-help]              0.04 0.04 -0.02 0.01 0.03 0.06  0.14     1930
#>                                 Tail_ESS Rhat
#> marg[1: Group counselling]          2362    1
#> marg[1: Individual counselling]     2103    1
#> marg[1: Self-help]                  2354    1
#> 
#> ---------------------------------------------------------------------- Study: 2 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[2: Group counselling]      0.13 0.08  0.02 0.07 0.11 0.17  0.31     2783
#> marg[2: Individual counselling] 0.09 0.05  0.02 0.05 0.08 0.11  0.21     2249
#> marg[2: Self-help]              0.05 0.05 -0.03 0.02 0.04 0.07  0.18     2190
#>                                 Tail_ESS Rhat
#> marg[2: Group counselling]          3016    1
#> marg[2: Individual counselling]     2457    1
#> marg[2: Self-help]                  2308    1
#> 
#> ---------------------------------------------------------------------- Study: 3 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[3: Group counselling]      0.16 0.09  0.03 0.10 0.15 0.21  0.37     1946
#> marg[3: Individual counselling] 0.11 0.04  0.04 0.08 0.11 0.14  0.21     1182
#> marg[3: Self-help]              0.06 0.06 -0.03 0.02 0.06 0.10  0.20     1838
#>                                 Tail_ESS Rhat
#> marg[3: Group counselling]          2373    1
#> marg[3: Individual counselling]     1656    1
#> marg[3: Self-help]                  2451    1
#> 
#> ---------------------------------------------------------------------- Study: 4 ---- 
#> 
#>                                 mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[4: Group counselling]      0.04 0.03 0.00 0.02 0.03 0.05  0.13     2825
#> marg[4: Individual counselling] 0.02 0.02 0.01 0.01 0.02 0.03  0.07     2553
#> marg[4: Self-help]              0.01 0.02 0.00 0.00 0.01 0.02  0.06     2064
#>                                 Tail_ESS Rhat
#> marg[4: Group counselling]          2624    1
#> marg[4: Individual counselling]     2771    1
#> marg[4: Self-help]                  2507    1
#> 
#> ---------------------------------------------------------------------- Study: 5 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[5: Group counselling]      0.16 0.09  0.03 0.10 0.15 0.21  0.37     1924
#> marg[5: Individual counselling] 0.11 0.04  0.04 0.08 0.11 0.14  0.21     1141
#> marg[5: Self-help]              0.06 0.06 -0.03 0.02 0.06 0.09  0.20     1842
#>                                 Tail_ESS Rhat
#> marg[5: Group counselling]          2246    1
#> marg[5: Individual counselling]     1825    1
#> marg[5: Self-help]                  2428    1
#> 
#> ---------------------------------------------------------------------- Study: 6 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[6: Group counselling]      0.07 0.06  0.01 0.03 0.06 0.10  0.23     3016
#> marg[6: Individual counselling] 0.05 0.03  0.01 0.02 0.04 0.06  0.12     2459
#> marg[6: Self-help]              0.03 0.03 -0.01 0.01 0.02 0.04  0.10     2388
#>                                 Tail_ESS Rhat
#> marg[6: Group counselling]          2569    1
#> marg[6: Individual counselling]     2428    1
#> marg[6: Self-help]                  2638    1
#> 
#> ---------------------------------------------------------------------- Study: 7 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[7: Group counselling]      0.09 0.06  0.01 0.05 0.08 0.12  0.24     2307
#> marg[7: Individual counselling] 0.06 0.03  0.02 0.04 0.05 0.07  0.13     1765
#> marg[7: Self-help]              0.03 0.03 -0.01 0.01 0.03 0.05  0.11     2124
#>                                 Tail_ESS Rhat
#> marg[7: Group counselling]          2549    1
#> marg[7: Individual counselling]     1952    1
#> marg[7: Self-help]                  2559    1
#> 
#> ---------------------------------------------------------------------- Study: 8 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[8: Group counselling]      0.12 0.08  0.01 0.06 0.10 0.16  0.30     2416
#> marg[8: Individual counselling] 0.08 0.04  0.02 0.05 0.07 0.10  0.18     2124
#> marg[8: Self-help]              0.04 0.04 -0.02 0.01 0.04 0.06  0.15     2238
#>                                 Tail_ESS Rhat
#> marg[8: Group counselling]          2670    1
#> marg[8: Individual counselling]     1952    1
#> marg[8: Self-help]                  2527    1
#> 
#> ---------------------------------------------------------------------- Study: 9 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[9: Group counselling]      0.19 0.10  0.03 0.12 0.18 0.25  0.42     2149
#> marg[9: Individual counselling] 0.13 0.05  0.04 0.09 0.13 0.17  0.25     1504
#> marg[9: Self-help]              0.08 0.07 -0.03 0.03 0.07 0.11  0.24     1849
#>                                 Tail_ESS Rhat
#> marg[9: Group counselling]          2404    1
#> marg[9: Individual counselling]     1868    1
#> marg[9: Self-help]                  2254    1
#> 
#> --------------------------------------------------------------------- Study: 10 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[10: Group counselling]      0.17 0.09  0.03 0.11 0.16 0.22  0.38     1945
#> marg[10: Individual counselling] 0.12 0.04  0.04 0.09 0.11 0.14  0.22     1191
#> marg[10: Self-help]              0.07 0.06 -0.03 0.03 0.06 0.10  0.21     1818
#>                                  Tail_ESS Rhat
#> marg[10: Group counselling]          2314    1
#> marg[10: Individual counselling]     1565    1
#> marg[10: Self-help]                  2413    1
#> 
#> --------------------------------------------------------------------- Study: 11 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[11: Group counselling]      0.05 0.04  0.01 0.03 0.05 0.07  0.15     2058
#> marg[11: Individual counselling] 0.03 0.02  0.01 0.02 0.03 0.04  0.07     1382
#> marg[11: Self-help]              0.02 0.02 -0.01 0.01 0.02 0.03  0.06     1858
#>                                  Tail_ESS Rhat
#> marg[11: Group counselling]          2223    1
#> marg[11: Individual counselling]     1753    1
#> marg[11: Self-help]                  2368    1
#> 
#> --------------------------------------------------------------------- Study: 12 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[12: Group counselling]      0.16 0.09  0.03 0.10 0.15 0.20  0.36     1954
#> marg[12: Individual counselling] 0.11 0.04  0.04 0.08 0.10 0.13  0.20     1168
#> marg[12: Self-help]              0.06 0.05 -0.02 0.02 0.05 0.09  0.19     1830
#>                                  Tail_ESS Rhat
#> marg[12: Group counselling]          2237    1
#> marg[12: Individual counselling]     1755    1
#> marg[12: Self-help]                  2383    1
#> 
#> --------------------------------------------------------------------- Study: 13 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[13: Group counselling]      0.12 0.08  0.02 0.06 0.10 0.16  0.31     2378
#> marg[13: Individual counselling] 0.08 0.04  0.02 0.05 0.07 0.10  0.17     1765
#> marg[13: Self-help]              0.04 0.04 -0.02 0.01 0.04 0.07  0.16     1982
#>                                  Tail_ESS Rhat
#> marg[13: Group counselling]          2560    1
#> marg[13: Individual counselling]     2326    1
#> marg[13: Self-help]                  2368    1
#> 
#> --------------------------------------------------------------------- Study: 14 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[14: Group counselling]      0.14 0.08  0.02 0.08 0.13 0.18  0.33     2016
#> marg[14: Individual counselling] 0.09 0.04  0.03 0.07 0.09 0.11  0.18     1281
#> marg[14: Self-help]              0.05 0.05 -0.02 0.02 0.05 0.08  0.17     1893
#>                                  Tail_ESS Rhat
#> marg[14: Group counselling]          2329    1
#> marg[14: Individual counselling]     1737    1
#> marg[14: Self-help]                  2241    1
#> 
#> --------------------------------------------------------------------- Study: 15 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[15: Group counselling]      0.11 0.07  0.01 0.06 0.10 0.15  0.29     2520
#> marg[15: Individual counselling] 0.08 0.05  0.01 0.04 0.07 0.11  0.19     2294
#> marg[15: Self-help]              0.04 0.05 -0.02 0.01 0.03 0.07  0.16     2297
#>                                  Tail_ESS Rhat
#> marg[15: Group counselling]          2464    1
#> marg[15: Individual counselling]     2769    1
#> marg[15: Self-help]                  2292    1
#> 
#> --------------------------------------------------------------------- Study: 16 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[16: Group counselling]      0.12 0.07  0.02 0.07 0.11 0.16  0.30     2129
#> marg[16: Individual counselling] 0.08 0.04  0.02 0.05 0.08 0.10  0.17     1600
#> marg[16: Self-help]              0.04 0.04 -0.02 0.02 0.04 0.07  0.15     1896
#>                                  Tail_ESS Rhat
#> marg[16: Group counselling]          2582    1
#> marg[16: Individual counselling]     2137    1
#> marg[16: Self-help]                  2460    1
#> 
#> --------------------------------------------------------------------- Study: 17 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[17: Group counselling]      0.14 0.08  0.02 0.09 0.13 0.18  0.33     1948
#> marg[17: Individual counselling] 0.09 0.04  0.03 0.07 0.09 0.12  0.18     1157
#> marg[17: Self-help]              0.05 0.05 -0.02 0.02 0.05 0.08  0.17     1830
#>                                  Tail_ESS Rhat
#> marg[17: Group counselling]          2276    1
#> marg[17: Individual counselling]     1657    1
#> marg[17: Self-help]                  2303    1
#> 
#> --------------------------------------------------------------------- Study: 18 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[18: Group counselling]      0.13 0.08  0.02 0.07 0.11 0.17  0.31     2075
#> marg[18: Individual counselling] 0.08 0.04  0.03 0.06 0.08 0.10  0.17     1362
#> marg[18: Self-help]              0.05 0.05 -0.02 0.02 0.04 0.07  0.16     1867
#>                                  Tail_ESS Rhat
#> marg[18: Group counselling]          2439    1
#> marg[18: Individual counselling]     1985    1
#> marg[18: Self-help]                  2401    1
#> 
#> --------------------------------------------------------------------- Study: 19 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[19: Group counselling]      0.19 0.09  0.03 0.12 0.18 0.24  0.40     1944
#> marg[19: Individual counselling] 0.13 0.05  0.05 0.10 0.13 0.16  0.24     1175
#> marg[19: Self-help]              0.07 0.07 -0.03 0.03 0.07 0.11  0.23     1836
#>                                  Tail_ESS Rhat
#> marg[19: Group counselling]          2424    1
#> marg[19: Individual counselling]     1564    1
#> marg[19: Self-help]                  2389    1
#> 
#> --------------------------------------------------------------------- Study: 20 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[20: Group counselling]      0.11 0.06  0.02 0.06 0.10 0.14  0.27     1972
#> marg[20: Individual counselling] 0.07 0.03  0.02 0.05 0.06 0.08  0.14     1196
#> marg[20: Self-help]              0.04 0.04 -0.01 0.01 0.03 0.06  0.13     1855
#>                                  Tail_ESS Rhat
#> marg[20: Group counselling]          2354    1
#> marg[20: Individual counselling]     1594    1
#> marg[20: Self-help]                  2324    1
#> 
#> --------------------------------------------------------------------- Study: 21 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[21: Group counselling]      0.22 0.10  0.04 0.15 0.22 0.29  0.43     2375
#> marg[21: Individual counselling] 0.17 0.06  0.05 0.12 0.17 0.21  0.29     1640
#> marg[21: Self-help]              0.09 0.08 -0.06 0.04 0.09 0.15  0.26     2108
#>                                  Tail_ESS Rhat
#> marg[21: Group counselling]          2673    1
#> marg[21: Individual counselling]     2167    1
#> marg[21: Self-help]                  2374    1
#> 
#> --------------------------------------------------------------------- Study: 22 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[22: Group counselling]      0.13 0.08  0.02 0.08 0.12 0.18  0.32     3794
#> marg[22: Individual counselling] 0.10 0.06  0.02 0.05 0.09 0.13  0.23     2497
#> marg[22: Self-help]              0.05 0.05 -0.03 0.02 0.04 0.07  0.17     2473
#>                                  Tail_ESS Rhat
#> marg[22: Group counselling]          2982    1
#> marg[22: Individual counselling]     2698    1
#> marg[22: Self-help]                  2459    1
#> 
#> --------------------------------------------------------------------- Study: 23 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[23: Group counselling]      0.14 0.09  0.02 0.08 0.13 0.19  0.34     2924
#> marg[23: Individual counselling] 0.10 0.06  0.02 0.06 0.09 0.13  0.23     2691
#> marg[23: Self-help]              0.06 0.06 -0.03 0.02 0.04 0.08  0.21     2197
#>                                  Tail_ESS Rhat
#> marg[23: Group counselling]          2522    1
#> marg[23: Individual counselling]     2613    1
#> marg[23: Self-help]                  2419    1
#> 
#> --------------------------------------------------------------------- Study: 24 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[24: Group counselling]      0.11 0.08  0.01 0.05 0.09 0.15  0.30     3287
#> marg[24: Individual counselling] 0.07 0.05  0.01 0.04 0.06 0.10  0.20     3106
#> marg[24: Self-help]              0.04 0.05 -0.02 0.01 0.03 0.06  0.17     2143
#>                                  Tail_ESS Rhat
#> marg[24: Group counselling]          2843    1
#> marg[24: Individual counselling]     2566    1
#> marg[24: Self-help]                  2198    1
#> 

# Since there are no covariates in the model, the marginal and conditional
# (log) odds ratios here coincide
marginal_effects(smk_fit_RE, mtype = "link")
#> ---------------------------------------------------------------------- Study: 1 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[1: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[1: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[1: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[1: Group counselling]          2375    1
#> marg[1: Individual counselling]     1549    1
#> marg[1: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 2 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[2: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[2: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[2: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[2: Group counselling]          2375    1
#> marg[2: Individual counselling]     1549    1
#> marg[2: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 3 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[3: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[3: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[3: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[3: Group counselling]          2375    1
#> marg[3: Individual counselling]     1549    1
#> marg[3: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 4 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[4: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[4: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[4: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[4: Group counselling]          2375    1
#> marg[4: Individual counselling]     1549    1
#> marg[4: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 5 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[5: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[5: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[5: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[5: Group counselling]          2375    1
#> marg[5: Individual counselling]     1549    1
#> marg[5: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 6 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[6: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[6: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[6: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[6: Group counselling]          2375    1
#> marg[6: Individual counselling]     1549    1
#> marg[6: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 7 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[7: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[7: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[7: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[7: Group counselling]          2375    1
#> marg[7: Individual counselling]     1549    1
#> marg[7: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 8 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[8: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[8: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[8: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[8: Group counselling]          2375    1
#> marg[8: Individual counselling]     1549    1
#> marg[8: Self-help]                  2353    1
#> 
#> ---------------------------------------------------------------------- Study: 9 ---- 
#> 
#>                                 mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[9: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[9: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[9: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                 Tail_ESS Rhat
#> marg[9: Group counselling]          2375    1
#> marg[9: Individual counselling]     1549    1
#> marg[9: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 10 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[10: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[10: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[10: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[10: Group counselling]          2375    1
#> marg[10: Individual counselling]     1549    1
#> marg[10: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 11 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[11: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[11: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[11: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[11: Group counselling]          2375    1
#> marg[11: Individual counselling]     1549    1
#> marg[11: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 12 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[12: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[12: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[12: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[12: Group counselling]          2375    1
#> marg[12: Individual counselling]     1549    1
#> marg[12: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 13 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[13: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[13: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[13: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[13: Group counselling]          2375    1
#> marg[13: Individual counselling]     1549    1
#> marg[13: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 14 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[14: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[14: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[14: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[14: Group counselling]          2375    1
#> marg[14: Individual counselling]     1549    1
#> marg[14: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 15 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[15: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[15: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[15: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[15: Group counselling]          2375    1
#> marg[15: Individual counselling]     1549    1
#> marg[15: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 16 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[16: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[16: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[16: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[16: Group counselling]          2375    1
#> marg[16: Individual counselling]     1549    1
#> marg[16: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 17 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[17: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[17: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[17: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[17: Group counselling]          2375    1
#> marg[17: Individual counselling]     1549    1
#> marg[17: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 18 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[18: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[18: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[18: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[18: Group counselling]          2375    1
#> marg[18: Individual counselling]     1549    1
#> marg[18: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 19 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[19: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[19: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[19: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[19: Group counselling]          2375    1
#> marg[19: Individual counselling]     1549    1
#> marg[19: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 20 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[20: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[20: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[20: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[20: Group counselling]          2375    1
#> marg[20: Individual counselling]     1549    1
#> marg[20: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 21 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[21: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[21: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[21: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[21: Group counselling]          2375    1
#> marg[21: Individual counselling]     1549    1
#> marg[21: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 22 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[22: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[22: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[22: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[22: Group counselling]          2375    1
#> marg[22: Individual counselling]     1549    1
#> marg[22: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 23 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[23: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[23: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[23: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[23: Group counselling]          2375    1
#> marg[23: Individual counselling]     1549    1
#> marg[23: Self-help]                  2353    1
#> 
#> --------------------------------------------------------------------- Study: 24 ---- 
#> 
#>                                  mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[24: Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> marg[24: Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> marg[24: Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                                  Tail_ESS Rhat
#> marg[24: Group counselling]          2375    1
#> marg[24: Individual counselling]     1549    1
#> marg[24: Self-help]                  2353    1
#> 
relative_effects(smk_fit_RE)
#>                           mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> d[Group counselling]      1.10 0.44  0.27 0.81 1.09 1.38  2.03     1924
#> d[Individual counselling] 0.84 0.25  0.37 0.67 0.84 1.00  1.35     1147
#> d[Self-help]              0.50 0.41 -0.30 0.24 0.49 0.76  1.30     1830
#>                           Tail_ESS Rhat
#> d[Group counselling]          2375    1
#> d[Individual counselling]     1549    1
#> d[Self-help]                  2353    1

# Marginal risk differences in a population with 67 observed events out of
# 566 individuals on No Intervention, corresponding to a Beta(67, 566 - 67)
# distribution on the baseline probability of response
(smk_rd_RE <- marginal_effects(smk_fit_RE,
                               baseline = distr(qbeta, 67, 566 - 67),
                               baseline_type = "response",
                               mtype = "difference"))
#>                              mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[Group counselling]      0.18 0.09  0.03 0.11 0.17 0.23  0.39     1922
#> marg[Individual counselling] 0.12 0.05  0.04 0.09 0.12 0.15  0.23     1190
#> marg[Self-help]              0.07 0.06 -0.03 0.03 0.06 0.10  0.21     1823
#>                              Tail_ESS Rhat
#> marg[Group counselling]          2263    1
#> marg[Individual counselling]     1625    1
#> marg[Self-help]                  2328    1
plot(smk_rd_RE)

# }

## Plaque psoriasis ML-NMR
# \donttest{
# Run plaque psoriasis ML-NMR example if not already available
if (!exists("pso_fit")) example("example_pso_mlnmr", run.donttest = TRUE)
# }
# \donttest{
# Population-average marginal probit differences in each study in the network
(pso_marg <- marginal_effects(pso_fit, mtype = "link"))
#> ---------------------------------------------------------------- Study: FIXTURE ---- 
#> 
#>                        mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[FIXTURE: ETN]     1.63 0.08 1.47 1.58 1.64 1.69  1.80     4908     3553
#> marg[FIXTURE: IXE_Q2W] 2.98 0.09 2.80 2.92 2.98 3.04  3.17     6443     3128
#> marg[FIXTURE: IXE_Q4W] 2.57 0.09 2.41 2.52 2.57 2.63  2.75     6464     3888
#> marg[FIXTURE: SEC_150] 2.18 0.11 1.96 2.11 2.18 2.26  2.41     5159     3524
#> marg[FIXTURE: SEC_300] 2.49 0.12 2.26 2.41 2.49 2.57  2.71     6028     3269
#>                        Rhat
#> marg[FIXTURE: ETN]        1
#> marg[FIXTURE: IXE_Q2W]    1
#> marg[FIXTURE: IXE_Q4W]    1
#> marg[FIXTURE: SEC_150]    1
#> marg[FIXTURE: SEC_300]    1
#> 
#> -------------------------------------------------------------- Study: UNCOVER-1 ---- 
#> 
#>                          mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[UNCOVER-1: ETN]     1.48 0.08 1.32 1.43 1.48 1.54  1.64     4769     3545
#> marg[UNCOVER-1: IXE_Q2W] 2.87 0.09 2.70 2.81 2.87 2.93  3.04     6614     3057
#> marg[UNCOVER-1: IXE_Q4W] 2.46 0.08 2.31 2.41 2.46 2.52  2.62     6891     3202
#> marg[UNCOVER-1: SEC_150] 2.07 0.12 1.83 1.99 2.07 2.16  2.30     5911     3547
#> marg[UNCOVER-1: SEC_300] 2.37 0.12 2.13 2.29 2.38 2.46  2.62     6765     3374
#>                          Rhat
#> marg[UNCOVER-1: ETN]        1
#> marg[UNCOVER-1: IXE_Q2W]    1
#> marg[UNCOVER-1: IXE_Q4W]    1
#> marg[UNCOVER-1: SEC_150]    1
#> marg[UNCOVER-1: SEC_300]    1
#> 
#> -------------------------------------------------------------- Study: UNCOVER-2 ---- 
#> 
#>                          mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[UNCOVER-2: ETN]     1.48 0.08 1.33 1.43 1.48 1.54  1.64     5272     3512
#> marg[UNCOVER-2: IXE_Q2W] 2.87 0.09 2.71 2.81 2.87 2.93  3.04     6991     2922
#> marg[UNCOVER-2: IXE_Q4W] 2.47 0.08 2.32 2.41 2.46 2.52  2.62     7495     3221
#> marg[UNCOVER-2: SEC_150] 2.08 0.12 1.83 2.00 2.07 2.15  2.30     5970     3604
#> marg[UNCOVER-2: SEC_300] 2.38 0.12 2.13 2.29 2.38 2.46  2.62     6890     3179
#>                          Rhat
#> marg[UNCOVER-2: ETN]        1
#> marg[UNCOVER-2: IXE_Q2W]    1
#> marg[UNCOVER-2: IXE_Q4W]    1
#> marg[UNCOVER-2: SEC_150]    1
#> marg[UNCOVER-2: SEC_300]    1
#> 
#> -------------------------------------------------------------- Study: UNCOVER-3 ---- 
#> 
#>                          mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[UNCOVER-3: ETN]     1.50 0.08 1.35 1.45 1.50 1.55  1.66     5310     3460
#> marg[UNCOVER-3: IXE_Q2W] 2.89 0.09 2.72 2.83 2.89 2.95  3.06     6992     2949
#> marg[UNCOVER-3: IXE_Q4W] 2.48 0.08 2.34 2.43 2.48 2.53  2.64     7557     3394
#> marg[UNCOVER-3: SEC_150] 2.09 0.11 1.86 2.02 2.09 2.17  2.32     5838     3537
#> marg[UNCOVER-3: SEC_300] 2.40 0.12 2.16 2.32 2.40 2.47  2.63     6786     3392
#>                          Rhat
#> marg[UNCOVER-3: ETN]        1
#> marg[UNCOVER-3: IXE_Q2W]    1
#> marg[UNCOVER-3: IXE_Q4W]    1
#> marg[UNCOVER-3: SEC_150]    1
#> marg[UNCOVER-3: SEC_300]    1
#> 
plot(pso_marg, ref_line = c(0, 1))


# Population-average marginal probit differences in a new target population,
# with means and SDs or proportions given by
new_agd_int <- data.frame(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1
)

# We need to add integration points to this data frame of new data
# We use the weighted mean correlation matrix computed from the IPD studies
new_agd_int <- add_integration(new_agd_int,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               cor = pso_net$int_cor,
                               n_int = 64)

# Population-average marginal probit differences of achieving PASI 75 in this
# target population, given a Normal(-1.75, 0.08^2) distribution on the
# baseline probit-probability of response on Placebo (at the reference levels
# of the covariates), are given by
(pso_marg_new <- marginal_effects(pso_fit,
                                  mtype = "link",
                                  newdata = new_agd_int,
                                  baseline = distr(qnorm, -1.75, 0.08)))
#> ------------------------------------------------------------------ Study: New 1 ---- 
#> 
#>                      mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> marg[New 1: ETN]     1.24 0.23 0.79 1.08 1.24 1.39  1.67     6555     3073    1
#> marg[New 1: IXE_Q2W] 2.85 0.21 2.43 2.71 2.86 3.00  3.26     7869     2961    1
#> marg[New 1: IXE_Q4W] 2.44 0.21 2.03 2.30 2.45 2.58  2.85     7732     2857    1
#> marg[New 1: SEC_150] 2.05 0.22 1.62 1.90 2.06 2.20  2.47     7488     3024    1
#> marg[New 1: SEC_300] 2.35 0.22 1.93 2.20 2.36 2.50  2.78     7713     2983    1
#> 
plot(pso_marg_new)

# }

## Progression free survival with newly-diagnosed multiple myeloma
# \donttest{
# Run newly-diagnosed multiple myeloma example if not already available
if (!exists("ndmm_fit")) example("example_ndmm", run.donttest = TRUE)
# }
# \donttest{
# We can produce a range of marginal effects from models with survival
# outcomes, specified with the mtype and type arguments. For example:

# Marginal survival probability difference at 5 years, all contrasts
marginal_effects(ndmm_fit, type = "survival", mtype = "difference",
                 times = 5, all_contrasts = TRUE)
#> -------------------------------------------------------------- Study: Attal2012 ---- 
#> 
#>                                  .time  mean   sd  2.5%   25%   50%   75% 97.5%
#> marg[Attal2012: Len vs. Pbo, 1]      5  0.19 0.02  0.16  0.18  0.19  0.20  0.22
#> marg[Attal2012: Thal vs. Pbo, 1]     5  0.04 0.03 -0.02  0.02  0.04  0.06  0.09
#> marg[Attal2012: Thal vs. Len, 1]     5 -0.15 0.03 -0.22 -0.17 -0.15 -0.13 -0.09
#>                                  Bulk_ESS Tail_ESS Rhat
#> marg[Attal2012: Len vs. Pbo, 1]      5315     2967    1
#> marg[Attal2012: Thal vs. Pbo, 1]     5100     3253    1
#> marg[Attal2012: Thal vs. Len, 1]     4954     3222    1
#> 
#> ------------------------------------------------------------ Study: Jackson2019 ---- 
#> 
#>                                    .time  mean   sd  2.5%   25%   50%   75%
#> marg[Jackson2019: Len vs. Pbo, 1]      5  0.19 0.02  0.16  0.18  0.19  0.21
#> marg[Jackson2019: Thal vs. Pbo, 1]     5  0.04 0.03 -0.02  0.02  0.04  0.06
#> marg[Jackson2019: Thal vs. Len, 1]     5 -0.16 0.03 -0.22 -0.18 -0.16 -0.13
#>                                    97.5% Bulk_ESS Tail_ESS Rhat
#> marg[Jackson2019: Len vs. Pbo, 1]   0.23     5313     3086    1
#> marg[Jackson2019: Thal vs. Pbo, 1]  0.10     5081     3252    1
#> marg[Jackson2019: Thal vs. Len, 1] -0.09     4954     3222    1
#> 
#> ----------------------------------------------------------- Study: McCarthy2012 ---- 
#> 
#>                                     .time  mean   sd  2.5%   25%   50%   75%
#> marg[McCarthy2012: Len vs. Pbo, 1]      5  0.20 0.02  0.16  0.18  0.19  0.21
#> marg[McCarthy2012: Thal vs. Pbo, 1]     5  0.04 0.03 -0.02  0.02  0.04  0.06
#> marg[McCarthy2012: Thal vs. Len, 1]     5 -0.15 0.03 -0.22 -0.18 -0.16 -0.13
#>                                     97.5% Bulk_ESS Tail_ESS Rhat
#> marg[McCarthy2012: Len vs. Pbo, 1]   0.23     5274     2996    1
#> marg[McCarthy2012: Thal vs. Pbo, 1]  0.10     5069     3281    1
#> marg[McCarthy2012: Thal vs. Len, 1] -0.09     4940     3181    1
#> 
#> ------------------------------------------------------------- Study: Morgan2012 ---- 
#> 
#>                                   .time  mean   sd  2.5%   25%   50%   75%
#> marg[Morgan2012: Len vs. Pbo, 1]      5  0.19 0.02  0.16  0.18  0.19  0.21
#> marg[Morgan2012: Thal vs. Pbo, 1]     5  0.04 0.03 -0.02  0.02  0.04  0.06
#> marg[Morgan2012: Thal vs. Len, 1]     5 -0.16 0.03 -0.22 -0.18 -0.16 -0.13
#>                                   97.5% Bulk_ESS Tail_ESS Rhat
#> marg[Morgan2012: Len vs. Pbo, 1]   0.23     5206     2933    1
#> marg[Morgan2012: Thal vs. Pbo, 1]  0.10     5058     3252    1
#> marg[Morgan2012: Thal vs. Len, 1] -0.09     4950     3224    1
#> 
#> ------------------------------------------------------------ Study: Palumbo2014 ---- 
#> 
#>                                    .time  mean   sd  2.5%   25%   50%   75%
#> marg[Palumbo2014: Len vs. Pbo, 1]      5  0.19 0.02  0.16  0.18  0.19  0.20
#> marg[Palumbo2014: Thal vs. Pbo, 1]     5  0.04 0.03 -0.02  0.02  0.04  0.06
#> marg[Palumbo2014: Thal vs. Len, 1]     5 -0.15 0.03 -0.21 -0.17 -0.15 -0.13
#>                                    97.5% Bulk_ESS Tail_ESS Rhat
#> marg[Palumbo2014: Len vs. Pbo, 1]   0.22     5318     2973    1
#> marg[Palumbo2014: Thal vs. Pbo, 1]  0.10     5063     3252    1
#> marg[Palumbo2014: Thal vs. Len, 1] -0.09     4984     3351    1
#> 

# Marginal difference in RMST up to 5 years
marginal_effects(ndmm_fit, type = "rmst", mtype = "difference", times = 5)
#> -------------------------------------------------------------- Study: Attal2012 ---- 
#> 
#>                          .time mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[Attal2012: Len, 1]      5 0.69 0.06  0.58 0.65 0.69 0.73  0.80     5271
#> marg[Attal2012: Thal, 1]     5 0.15 0.12 -0.08 0.07 0.15 0.23  0.37     5096
#>                          Tail_ESS Rhat
#> marg[Attal2012: Len, 1]      2981    1
#> marg[Attal2012: Thal, 1]     3281    1
#> 
#> ------------------------------------------------------------ Study: Jackson2019 ---- 
#> 
#>                            .time mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[Jackson2019: Len, 1]      5 0.74 0.06  0.62 0.70 0.74 0.78  0.86     5077
#> marg[Jackson2019: Thal, 1]     5 0.16 0.13 -0.09 0.08 0.16 0.25  0.40     5060
#>                            Tail_ESS Rhat
#> marg[Jackson2019: Len, 1]      2663    1
#> marg[Jackson2019: Thal, 1]     3252    1
#> 
#> ----------------------------------------------------------- Study: McCarthy2012 ---- 
#> 
#>                             .time mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[McCarthy2012: Len, 1]      5 0.64 0.05  0.53 0.60 0.64 0.67  0.75     4960
#> marg[McCarthy2012: Thal, 1]     5 0.14 0.11 -0.08 0.07 0.14 0.22  0.35     5078
#>                             Tail_ESS Rhat
#> marg[McCarthy2012: Len, 1]      2863    1
#> marg[McCarthy2012: Thal, 1]     3252    1
#> 
#> ------------------------------------------------------------- Study: Morgan2012 ---- 
#> 
#>                           .time mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[Morgan2012: Len, 1]      5 0.76 0.06  0.65 0.72 0.76 0.80  0.88     5319
#> marg[Morgan2012: Thal, 1]     5 0.17 0.13 -0.09 0.08 0.17 0.26  0.42     5075
#>                           Tail_ESS Rhat
#> marg[Morgan2012: Len, 1]      2814    1
#> marg[Morgan2012: Thal, 1]     3252    1
#> 
#> ------------------------------------------------------------ Study: Palumbo2014 ---- 
#> 
#>                            .time mean   sd  2.5%  25%  50%  75% 97.5% Bulk_ESS
#> marg[Palumbo2014: Len, 1]      5 0.75 0.07  0.63 0.71 0.75 0.80  0.88     5090
#> marg[Palumbo2014: Thal, 1]     5 0.16 0.13 -0.09 0.08 0.16 0.25  0.40     5078
#>                            Tail_ESS Rhat
#> marg[Palumbo2014: Len, 1]      2935    1
#> marg[Palumbo2014: Thal, 1]     3281    1
#> 

# Marginal median survival time ratios
marginal_effects(ndmm_fit, type = "median", mtype = "ratio")
#> -------------------------------------------------------------- Study: Attal2012 ---- 
#> 
#>                       mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[Attal2012: Len]  1.52 0.06 1.40 1.48 1.51 1.55  1.64     5118     3195
#> marg[Attal2012: Thal] 1.09 0.07 0.96 1.04 1.09 1.14  1.24     5134     3311
#>                       Rhat
#> marg[Attal2012: Len]     1
#> marg[Attal2012: Thal]    1
#> 
#> ------------------------------------------------------------ Study: Jackson2019 ---- 
#> 
#>                         mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[Jackson2019: Len]  1.78 0.09 1.62 1.72 1.78 1.84  1.96     4849     2974
#> marg[Jackson2019: Thal] 1.13 0.11 0.94 1.06 1.13 1.20  1.35     5069     3241
#>                         Rhat
#> marg[Jackson2019: Len]     1
#> marg[Jackson2019: Thal]    1
#> 
#> ----------------------------------------------------------- Study: McCarthy2012 ---- 
#> 
#>                          mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[McCarthy2012: Len]  1.52 0.06 1.41 1.48 1.52 1.56  1.65     4822     2547
#> marg[McCarthy2012: Thal] 1.09 0.07 0.95 1.04 1.09 1.14  1.24     5025     3227
#>                          Rhat
#> marg[McCarthy2012: Len]     1
#> marg[McCarthy2012: Thal]    1
#> 
#> ------------------------------------------------------------- Study: Morgan2012 ---- 
#> 
#>                        mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[Morgan2012: Len]  1.85 0.10 1.66 1.78 1.84 1.91  2.06     5175     3289
#> marg[Morgan2012: Thal] 1.14 0.11 0.93 1.06 1.14 1.22  1.37     5049     3252
#>                        Rhat
#> marg[Morgan2012: Len]     1
#> marg[Morgan2012: Thal]    1
#> 
#> ------------------------------------------------------------ Study: Palumbo2014 ---- 
#> 
#>                         mean  sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS
#> marg[Palumbo2014: Len]  1.71 0.1 1.52 1.64 1.70 1.77  1.91     4727     3094
#> marg[Palumbo2014: Thal] 1.12 0.1 0.94 1.05 1.12 1.18  1.32     5010     3281
#>                         Rhat
#> marg[Palumbo2014: Len]     1
#> marg[Palumbo2014: Thal]    1
#> 

# Marginal log hazard ratios
# With no covariates in the model, these are constant over time and study
# populations, and are equal to the log hazard ratios from relative_effects()
plot(marginal_effects(ndmm_fit, type = "hazard", mtype = "link"),
     # The hazard is infinite at t=0 in some studies, giving undefined logHRs at t=0
     na.rm = TRUE)


# The NDMM vignette demonstrates the production of time-varying marginal
# hazard ratios from a ML-NMR model that includes covariates, see
# `vignette("example_ndmm")`

# Marginal survival difference over time
plot(marginal_effects(ndmm_fit, type = "survival", mtype = "difference"))

# }
```
