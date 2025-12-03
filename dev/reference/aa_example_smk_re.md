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
#> d[Group counselling]          1.11    0.01 0.45     0.26     0.82     1.10
#> d[Individual counselling]     0.85    0.01 0.24     0.39     0.69     0.84
#> d[Self-help]                  0.48    0.01 0.41    -0.30     0.22     0.48
#> lp__                      -5767.79    0.19 6.42 -5781.22 -5771.84 -5767.54
#> tau                           0.84    0.01 0.19     0.55     0.71     0.82
#>                                75%    97.5% n_eff Rhat
#> d[Group counselling]          1.39     2.05  1472    1
#> d[Individual counselling]     1.01     1.36   969    1
#> d[Self-help]                  0.74     1.32  1472    1
#> lp__                      -5763.28 -5756.12  1126    1
#> tau                           0.95     1.27  1089    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Dec  3 16:25:35 2025.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
