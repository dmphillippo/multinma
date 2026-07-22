# Example newly-diagnosed multiple myeloma

Calling `example("example_ndmm")` will run a proportional hazards
Weibull NMA model on the newly-diagnosed multiple myeloma data, using
the code in the Examples section below. The resulting `stan_nma` object
`ndmm_fit` will then be available in the global environment.

## Examples

``` r
# Set up newly-diagnosed multiple myeloma network

head(ndmm_ipd)
#>          study trt       studyf trtf      age iss_stage3 response_cr_vgpr male
#> 1 McCarthy2012 Pbo McCarthy2012  Pbo 50.81625          0                1    0
#> 2 McCarthy2012 Pbo McCarthy2012  Pbo 62.18165          0                0    0
#> 3 McCarthy2012 Pbo McCarthy2012  Pbo 51.53762          1                1    1
#> 4 McCarthy2012 Pbo McCarthy2012  Pbo 46.74128          0                1    1
#> 5 McCarthy2012 Pbo McCarthy2012  Pbo 62.62561          0                1    1
#> 6 McCarthy2012 Pbo McCarthy2012  Pbo 49.24520          1                1    0
#>   eventtime status
#> 1 31.106516      1
#> 2  3.299623      0
#> 3 57.400000      0
#> 4 57.400000      0
#> 5 57.400000      0
#> 6 30.714460      0
head(ndmm_agd)
#>        study     studyf trt trtf eventtime status
#> 1 Morgan2012 Morgan2012 Pbo  Pbo  18.72575      1
#> 2 Morgan2012 Morgan2012 Pbo  Pbo  63.36000      0
#> 3 Morgan2012 Morgan2012 Pbo  Pbo  34.35726      1
#> 4 Morgan2012 Morgan2012 Pbo  Pbo  10.77826      1
#> 5 Morgan2012 Morgan2012 Pbo  Pbo  63.36000      0
#> 6 Morgan2012 Morgan2012 Pbo  Pbo  14.52966      1

ndmm_net <- combine_network(
  set_ipd(ndmm_ipd,
          study, trt,
          Surv = Surv(eventtime / 12, status)),
  set_agd_surv(ndmm_agd,
               study, trt,
               Surv = Surv(eventtime / 12, status),
               covariates = ndmm_agd_covs))
# \donttest{
# Fit Weibull (PH) model
ndmm_fit <- nma(ndmm_net, 
                likelihood = "weibull",
                prior_intercept = normal(scale = 100),
                prior_trt = normal(scale = 10),
                prior_aux = half_normal(scale = 10))
#> Note: Setting "Pbo" as the network reference treatment.

ndmm_fit
#> A fixed effects NMA with a weibull likelihood (log link).
#> Inference for Stan model: survival_param.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                         mean se_mean   sd     2.5%      25%      50%      75%
#> d[Len]                 -0.54    0.00 0.04    -0.63    -0.57    -0.54    -0.51
#> d[Thal]                -0.11    0.00 0.09    -0.28    -0.17    -0.11    -0.05
#> lp__                -6230.03    0.06 2.43 -6235.71 -6231.45 -6229.77 -6228.25
#> shape[Attal2012]        1.30    0.00 0.06     1.18     1.25     1.30     1.34
#> shape[Jackson2019]      0.93    0.00 0.02     0.89     0.92     0.93     0.95
#> shape[McCarthy2012]     1.29    0.00 0.07     1.16     1.25     1.29     1.34
#> shape[Morgan2012]       0.88    0.00 0.03     0.82     0.86     0.88     0.90
#> shape[Palumbo2014]      1.02    0.00 0.07     0.88     0.97     1.01     1.06
#>                        97.5% n_eff Rhat
#> d[Len]                 -0.45  5305    1
#> d[Thal]                 0.06  5012    1
#> lp__                -6226.20  1820    1
#> shape[Attal2012]        1.42  4170    1
#> shape[Jackson2019]      0.98  4719    1
#> shape[McCarthy2012]     1.43  4930    1
#> shape[Morgan2012]       0.94  4713    1
#> shape[Palumbo2014]      1.16  4706    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Jul 22 14:44:36 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
