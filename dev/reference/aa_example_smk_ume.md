# Example smoking UME NMA

Calling `example("example_smk_ume")` will run an unrelated mean effects
(inconsistency) NMA model with the smoking cessation data, using the
code in the Examples section below. The resulting `stan_nma` object
`smk_fit_RE_UME` will then be available in the global environment.

## Details

Smoking UME NMA for use in examples.

## Examples

``` r
# Set up network of smoking cessation data
head(smoking)
#>   studyn trtn                   trtc  r   n
#> 1      1    1        No intervention  9 140
#> 2      1    3 Individual counselling 23 140
#> 3      1    4      Group counselling 10 138
#> 4      2    2              Self-help 11  78
#> 5      2    3 Individual counselling 12  85
#> 6      2    4      Group counselling 29 170

smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Print details
smk_net
#> A network with 24 AgD studies (arm-based).
#> 
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study Treatment arms                                                 
#>  1     3: No intervention | Group counselling | Individual counselling
#>  2     3: Group counselling | Individual counselling | Self-help      
#>  3     2: No intervention | Individual counselling                    
#>  4     2: No intervention | Individual counselling                    
#>  5     2: No intervention | Individual counselling                    
#>  6     2: No intervention | Individual counselling                    
#>  7     2: No intervention | Individual counselling                    
#>  8     2: No intervention | Individual counselling                    
#>  9     2: No intervention | Individual counselling                    
#>  10    2: No intervention | Self-help                                 
#>  ... plus 14 more studies
#> 
#>  Outcome type: count
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 4
#> Total number of studies: 24
#> Reference treatment is: No intervention
#> Network is connected

# \donttest{
# Fitting an unrelated mean effects (inconsistency) model
smk_fit_RE_UME <- nma(smk_net, 
                      consistency = "ume",
                      trt_effects = "random",
                      prior_intercept = normal(scale = 100),
                      prior_trt = normal(scale = 100),
                      prior_het = normal(scale = 5))

smk_fit_RE_UME
#> A random effects NMA with a binomial likelihood (logit link).
#> An inconsistency model ('ume') was fitted.
#> Inference for Stan model: binomial_1par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                                                     mean se_mean   sd     2.5%
#> d[Group counselling vs. No intervention]            1.16    0.02 0.80    -0.33
#> d[Individual counselling vs. No intervention]       0.91    0.01 0.29     0.34
#> d[Self-help vs. No intervention]                    0.34    0.01 0.61    -0.89
#> d[Individual counselling vs. Group counselling]    -0.28    0.01 0.63    -1.50
#> d[Self-help vs. Group counselling]                 -0.63    0.01 0.70    -2.01
#> d[Self-help vs. Individual counselling]             0.17    0.02 1.08    -1.92
#> lp__                                            -5765.51    0.19 6.41 -5778.95
#> tau                                                 0.94    0.01 0.23     0.59
#>                                                      25%      50%      75%
#> d[Group counselling vs. No intervention]            0.63     1.12     1.66
#> d[Individual counselling vs. No intervention]       0.72     0.90     1.09
#> d[Self-help vs. No intervention]                   -0.04     0.34     0.71
#> d[Individual counselling vs. Group counselling]    -0.69    -0.29     0.12
#> d[Self-help vs. Group counselling]                 -1.09    -0.63    -0.16
#> d[Self-help vs. Individual counselling]            -0.53     0.17     0.87
#> lp__                                            -5769.62 -5765.17 -5760.90
#> tau                                                 0.78     0.91     1.07
#>                                                    97.5% n_eff Rhat
#> d[Group counselling vs. No intervention]            2.80  2333 1.00
#> d[Individual counselling vs. No intervention]       1.48   738 1.00
#> d[Self-help vs. No intervention]                    1.54  1653 1.00
#> d[Individual counselling vs. Group counselling]     1.00  2183 1.00
#> d[Self-help vs. Group counselling]                  0.76  2304 1.00
#> d[Self-help vs. Individual counselling]             2.27  2968 1.00
#> lp__                                            -5754.05  1148 1.00
#> tau                                                 1.48   895 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 10 14:35:41 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
