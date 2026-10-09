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
#> d[Group counselling vs. No intervention]            1.13    0.01 0.79    -0.32
#> d[Individual counselling vs. No intervention]       0.91    0.01 0.28     0.39
#> d[Self-help vs. No intervention]                    0.35    0.01 0.59    -0.81
#> d[Individual counselling vs. Group counselling]    -0.29    0.01 0.60    -1.44
#> d[Self-help vs. Group counselling]                 -0.60    0.01 0.72    -2.01
#> d[Self-help vs. Individual counselling]             0.16    0.02 1.08    -2.06
#> lp__                                            -5765.24    0.19 6.17 -5778.26
#> tau                                                 0.93    0.01 0.22     0.59
#>                                                      25%      50%      75%
#> d[Group counselling vs. No intervention]            0.58     1.10     1.63
#> d[Individual counselling vs. No intervention]       0.73     0.91     1.09
#> d[Self-help vs. No intervention]                   -0.03     0.36     0.73
#> d[Individual counselling vs. Group counselling]    -0.68    -0.30     0.08
#> d[Self-help vs. Group counselling]                 -1.07    -0.58    -0.13
#> d[Self-help vs. Individual counselling]            -0.52     0.17     0.86
#> lp__                                            -5769.20 -5764.93 -5760.85
#> tau                                                 0.77     0.90     1.06
#>                                                    97.5% n_eff Rhat
#> d[Group counselling vs. No intervention]            2.78  2930    1
#> d[Individual counselling vs. No intervention]       1.50  1265    1
#> d[Self-help vs. No intervention]                    1.52  1993    1
#> d[Individual counselling vs. Group counselling]     0.93  2738    1
#> d[Self-help vs. Group counselling]                  0.79  2921    1
#> d[Self-help vs. Individual counselling]             2.31  3754    1
#> lp__                                            -5754.10  1031    1
#> tau                                                 1.43  1335    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Oct  9 09:15:34 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
