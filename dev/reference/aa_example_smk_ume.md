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
#> d[Group counselling vs. No intervention]            1.13    0.02 0.81    -0.40
#> d[Individual counselling vs. No intervention]       0.89    0.01 0.27     0.37
#> d[Self-help vs. No intervention]                    0.36    0.01 0.59    -0.81
#> d[Individual counselling vs. Group counselling]    -0.29    0.01 0.62    -1.51
#> d[Self-help vs. Group counselling]                 -0.62    0.01 0.71    -2.01
#> d[Self-help vs. Individual counselling]             0.16    0.02 1.07    -1.94
#> lp__                                            -5765.44    0.22 6.46 -5778.94
#> tau                                                 0.93    0.01 0.23     0.59
#>                                                      25%      50%      75%
#> d[Group counselling vs. No intervention]            0.60     1.10     1.63
#> d[Individual counselling vs. No intervention]       0.72     0.87     1.05
#> d[Self-help vs. No intervention]                    0.00     0.34     0.72
#> d[Individual counselling vs. Group counselling]    -0.68    -0.31     0.10
#> d[Self-help vs. Group counselling]                 -1.08    -0.61    -0.16
#> d[Self-help vs. Individual counselling]            -0.55     0.17     0.85
#> lp__                                            -5769.52 -5765.31 -5760.85
#> tau                                                 0.77     0.90     1.06
#>                                                    97.5% n_eff Rhat
#> d[Group counselling vs. No intervention]            2.80  2493    1
#> d[Individual counselling vs. No intervention]       1.45  1136    1
#> d[Self-help vs. No intervention]                    1.52  2415    1
#> d[Individual counselling vs. Group counselling]     0.96  2578    1
#> d[Self-help vs. Group counselling]                  0.79  2486    1
#> d[Self-help vs. Individual counselling]             2.24  3700    1
#> lp__                                            -5753.63   894    1
#> tau                                                 1.48  1129    1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Feb 24 14:08:53 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
