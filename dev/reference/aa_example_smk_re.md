# Example smoking RE NMA

Calling `example("example_smk_re")` will run a random effects NMA model
with the smoking cessation data, using the code in the Examples section
below. The resulting `stan_nma` object `smk_fit_RE` will then be
available in the global environment.

## Details

Smoking RE NMA for use in examples.

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
# Fitting a random effects model
smk_fit_RE <- nma(smk_net, 
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = normal(scale = 5))

smk_fit_RE
#> A random effects NMA with a binomial likelihood (logit link).
#> Inference for Stan model: binomial_1par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                               mean se_mean   sd     2.5%      25%      50%
#> d[Group counselling]          1.10    0.01 0.44     0.27     0.81     1.09
#> d[Individual counselling]     0.84    0.01 0.25     0.37     0.67     0.84
#> d[Self-help]                  0.50    0.01 0.41    -0.30     0.24     0.49
#> lp__                      -5768.43    0.21 6.44 -5781.69 -5772.60 -5768.21
#> tau                           0.84    0.01 0.19     0.54     0.71     0.81
#>                                75%    97.5% n_eff Rhat
#> d[Group counselling]          1.38     2.03  1879 1.00
#> d[Individual counselling]     1.00     1.35  1146 1.00
#> d[Self-help]                  0.76     1.30  1815 1.00
#> lp__                      -5763.92 -5756.55   972 1.00
#> tau                           0.94     1.27  1134 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Oct  8 10:45:45 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
