# Example smoking FE NMA

Calling `example("example_smk_fe")` will run a fixed effects NMA model
with the smoking cessation data, using the code in the Examples section
below. The resulting `stan_nma` object `smk_fit_FE` will then be
available in the global environment.

## Details

Smoking FE NMA for use in examples.

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
# Fitting a fixed effect model
smk_fit_FE <- nma(smk_net,
                  trt_effects = "fixed",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100))

smk_fit_FE
#> A fixed effects NMA with a binomial likelihood (logit link).
#> Inference for Stan model: binomial_1par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                               mean se_mean   sd     2.5%      25%      50%
#> d[Group counselling]          0.84    0.00 0.18     0.49     0.72     0.84
#> d[Individual counselling]     0.76    0.00 0.06     0.65     0.72     0.76
#> d[Self-help]                  0.22    0.00 0.13    -0.03     0.14     0.23
#> lp__                      -5859.40    0.09 3.64 -5867.59 -5861.60 -5859.08
#>                                75%    97.5% n_eff Rhat
#> d[Group counselling]          0.96     1.18  2302    1
#> d[Individual counselling]     0.80     0.88  1553    1
#> d[Self-help]                  0.31     0.47  2704    1
#> lp__                      -5856.84 -5853.18  1589    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Oct  9 11:09:32 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
