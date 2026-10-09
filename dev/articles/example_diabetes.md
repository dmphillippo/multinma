# Example: Diabetes

``` r

library(multinma)
options(mc.cores = parallel::detectCores())
```

    #> For execution on a local, multicore CPU with excess RAM we recommend calling
    #> options(mc.cores = parallel::detectCores())
    #> 
    #> Attaching package: 'multinma'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     dgamma, pgamma, qgamma

This vignette describes the analysis of data on the number of new cases
of diabetes in 22 trials of 6 antihypertensive drugs ([Elliott and Meyer
2007](#ref-Elliott2007); [Dias et al. 2011](#ref-TSD2)). The data are
available in this package as `diabetes`:

``` r

head(diabetes)
#>   studyn studyc trtn         trtc   r    n time
#> 1      1  MRC-E    1     Diuretic  43 1081  5.8
#> 2      1  MRC-E    2      Placebo  34 2213  5.8
#> 3      1  MRC-E    3 Beta Blocker  37 1102  5.8
#> 4      2   EWPH    1     Diuretic  29  416  4.7
#> 5      2   EWPH    2      Placebo  20  424  4.7
#> 6      3   SHEP    1     Diuretic 140 1631  3.0
```

## Setting up the network

We begin by setting up the network. We have arm-level count data giving
the number of new cases of diabetes (`r`) out of the total (`n`) in each
arm, so we use the function
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md).
For computational efficiency, we let “Beta Blocker” be set as the
network reference treatment by default. Elliott and Meyer
([2007](#ref-Elliott2007)) and Dias et al. ([2011](#ref-TSD2)) use
“Diuretic” as the reference, but it is a simple matter to transform the
results after fitting the NMA model.[^1]

``` r

db_net <- set_agd_arm(diabetes, 
                      study = studyc,
                      trt = trtc,
                      r = r, 
                      n = n)
db_net
#> A network with 22 AgD studies (arm-based).
#> 
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study  Treatment arms                       
#>  AASK   3: Beta Blocker | ACE Inhibitor | CCB
#>  ALLHAT 3: ACE Inhibitor | CCB | Diuretic    
#>  ALPINE 2: ARB | Diuretic                    
#>  ANBP-2 2: ACE Inhibitor | Diuretic          
#>  ASCOT  2: Beta Blocker | CCB                
#>  CAPPP  2: Beta Blocker | ACE Inhibitor      
#>  CHARM  2: ARB | Placebo                     
#>  DREAM  2: ACE Inhibitor | Placebo           
#>  EWPH   2: Diuretic | Placebo                
#>  FEVER  2: CCB | Placebo                     
#>  ... plus 12 more studies
#> 
#>  Outcome type: count
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 6
#> Total number of studies: 22
#> Reference treatment is: Beta Blocker
#> Network is connected
```

We also have details of length of follow-up in years in each trial
(`time`), which we will use as an offset with a cloglog link function to
model the data as rates. We do not have to specify this in the function
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md):
any additional columns in the data (e.g. offsets or covariates, here the
column `time`) will automatically be made available in the network.

Plot the network structure.

``` r

plot(db_net, weight_edges = TRUE, weight_nodes = TRUE)
```

![](example_diabetes_files/figure-html/diabetes_network_plot-1.png)

## Meta-analysis models

We fit both fixed effect (FE) and random effects (RE) models.

### Fixed effect meta-analysis

First, we fit a fixed effect model using the
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
function with `trt_effects = "fixed"`. We use \mathrm{N}(0, 100^2) prior
distributions for the treatment effects d_k and study-specific
intercepts \mu_j. We can examine the range of parameter values implied
by these prior distributions with the
[`summary()`](https://rdrr.io/r/base/summary.html) method:

``` r

summary(normal(scale = 100))
#> A Normal prior distribution: location = 0, scale = 100.
#> 50% of the prior density lies between -67.45 and 67.45.
#> 95% of the prior density lies between -196 and 196.
```

The model is fitted using the
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
function. We specify that a cloglog link will be used with
`link = "cloglog"` (the Binomial likelihood is the default for these
data), and specify the log follow-up time offset using the regression
formula `regression = ~offset(log(time))`.

``` r

db_fit_FE <- nma(db_net, 
                 trt_effects = "fixed",
                 link = "cloglog",
                 regression = ~offset(log(time)),
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100))
#> Note: No treatment classes specified in network, any interactions in `regression` formula will be separate (independent) for each treatment.
#> Use set_*() argument `trt_class` and nma() argument `class_interactions` to change this.
#> Note: Setting "Beta Blocker" as the network reference treatment.
```

Basic parameter summaries are given by the
[`print()`](https://rdrr.io/r/base/print.html) method:

``` r

db_fit_FE
#> A fixed effects NMA with a binomial likelihood (cloglog link).
#> Regression model: ~offset(log(time)).
#> Inference for Stan model: binomial_1par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                       mean se_mean   sd      2.5%       25%       50%       75%     97.5%
#> d[ACE Inhibitor]     -0.30    0.00 0.05     -0.39     -0.34     -0.30     -0.27     -0.21
#> d[ARB]               -0.40    0.00 0.05     -0.49     -0.43     -0.40     -0.37     -0.31
#> d[CCB]               -0.20    0.00 0.03     -0.26     -0.22     -0.20     -0.18     -0.14
#> d[Diuretic]           0.05    0.00 0.06     -0.05      0.02      0.06      0.09      0.16
#> d[Placebo]           -0.19    0.00 0.05     -0.29     -0.23     -0.19     -0.16     -0.09
#> lp__             -37970.30    0.09 3.64 -37978.50 -37972.56 -37970.06 -37967.69 -37963.89
#>                  n_eff Rhat
#> d[ACE Inhibitor]  1657    1
#> d[ARB]            2483    1
#> d[CCB]            1915    1
#> d[Diuretic]       1868    1
#> d[Placebo]        1482    1
#> lp__              1824    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Oct  9 09:35:43 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

By default, summaries of the study-specific intercepts \mu_j are hidden,
but could be examined by changing the `pars` argument:

``` r

# Not run
print(db_fit_FE, pars = c("d", "mu"))
```

The prior and posterior distributions can be compared visually using the
[`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/dev/reference/plot_prior_posterior.md)
function:

``` r

plot_prior_posterior(db_fit_FE)
```

![](example_diabetes_files/figure-html/db_FE_pp_plot-1.png)

### Random effects meta-analysis

We now fit a random effects model using the
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
function with `trt_effects = "random"`. Again, we use \mathrm{N}(0,
100^2) prior distributions for the treatment effects d_k and
study-specific intercepts \mu_j, and we additionally use a
\textrm{half-N}(5^2) prior for the heterogeneity standard deviation
\tau. We can examine the range of parameter values implied by these
prior distributions with the
[`summary()`](https://rdrr.io/r/base/summary.html) method:

``` r

summary(normal(scale = 100))
#> A Normal prior distribution: location = 0, scale = 100.
#> 50% of the prior density lies between -67.45 and 67.45.
#> 95% of the prior density lies between -196 and 196.
summary(half_normal(scale = 5))
#> A half-Normal prior distribution: scale = 5.
#> 50% of the prior density lies between 0 and 3.37.
#> 95% of the prior density lies between 0 and 9.8.
```

Fitting the RE model

``` r

db_fit_RE <- nma(db_net, 
                 trt_effects = "random",
                 link = "cloglog",
                 regression = ~offset(log(time)),
                 prior_intercept = normal(scale = 10),
                 prior_trt = normal(scale = 10),
                 prior_het = half_normal(scale = 5),
                 init_r = 0.5)
#> Note: No treatment classes specified in network, any interactions in `regression` formula will be separate (independent) for each treatment.
#> Use set_*() argument `trt_class` and nma() argument `class_interactions` to change this.
#> Note: Setting "Beta Blocker" as the network reference treatment.
```

Basic parameter summaries are given by the
[`print()`](https://rdrr.io/r/base/print.html) method:

``` r

db_fit_RE
#> A random effects NMA with a binomial likelihood (cloglog link).
#> Regression model: ~offset(log(time)).
#> Inference for Stan model: binomial_1par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                       mean se_mean   sd      2.5%       25%       50%       75%     97.5%
#> d[ACE Inhibitor]     -0.33    0.00 0.08     -0.49     -0.38     -0.33     -0.28     -0.19
#> d[ARB]               -0.40    0.00 0.10     -0.60     -0.46     -0.40     -0.34     -0.22
#> d[CCB]               -0.17    0.00 0.06     -0.30     -0.21     -0.17     -0.13     -0.04
#> d[Diuretic]           0.07    0.00 0.09     -0.10      0.01      0.07      0.13      0.24
#> d[Placebo]           -0.22    0.00 0.08     -0.39     -0.27     -0.21     -0.16     -0.06
#> lp__             -37980.54    0.22 6.96 -37994.94 -37985.04 -37980.25 -37975.70 -37968.10
#> tau                   0.13    0.00 0.04      0.05      0.10      0.13      0.16      0.23
#>                  n_eff Rhat
#> d[ACE Inhibitor]  1875    1
#> d[ARB]            2142    1
#> d[CCB]            2161    1
#> d[Diuretic]       2100    1
#> d[Placebo]        1759    1
#> lp__               992    1
#> tau                999    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Oct  9 09:35:58 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

By default, summaries of the study-specific intercepts \mu_j and
study-specific relative effects \delta\_{jk} are hidden, but could be
examined by changing the `pars` argument:

``` r

# Not run
print(db_fit_RE, pars = c("d", "mu", "delta"))
```

The prior and posterior distributions can be compared visually using the
[`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/dev/reference/plot_prior_posterior.md)
function:

``` r

plot_prior_posterior(db_fit_RE, prior = c("trt", "het"))
```

![](example_diabetes_files/figure-html/db_RE_pp_plot-1.png)

### Model comparison

Model fit can be checked using the
[`dic()`](https://dmphillippo.github.io/multinma/dev/reference/dic.md)
function:

``` r

(dic_FE <- dic(db_fit_FE))
#> Residual deviance: 78.2 (on 48 data points)
#>                pD: 27
#>               DIC: 105.2
```

``` r

(dic_RE <- dic(db_fit_RE))
#> Residual deviance: 53.3 (on 48 data points)
#>                pD: 37.9
#>               DIC: 91.2
```

The FE model is a very poor fit to the data, with a residual deviance
much higher than the number of data points. The RE model fits the data
better, and has a much lower DIC; we prefer the RE model.

We can also examine the residual deviance contributions with the
corresponding [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method.

``` r

plot(dic_FE)
```

![](example_diabetes_files/figure-html/db_FE_resdev_plot-1.png)

``` r

plot(dic_RE)
```

![](example_diabetes_files/figure-html/db_RE_resdev_plot-1.png)

## Further results

For comparison with Elliott and Meyer ([2007](#ref-Elliott2007)) and
Dias et al. ([2011](#ref-TSD2)), we can produce relative effects against
“Diuretic” using the
[`relative_effects()`](https://dmphillippo.github.io/multinma/dev/reference/relative_effects.md)
function with `trt_ref = "Diuretic"`:

``` r

(db_releff_FE <- relative_effects(db_fit_FE, trt_ref = "Diuretic"))
#>                   mean   sd  2.5%   25%   50%   75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d[Beta Blocker]  -0.05 0.06 -0.16 -0.09 -0.06 -0.02  0.05     1895     2372    1
#> d[ACE Inhibitor] -0.36 0.05 -0.46 -0.39 -0.36 -0.32 -0.26     5041     3299    1
#> d[ARB]           -0.45 0.06 -0.58 -0.50 -0.45 -0.41 -0.33     4056     3578    1
#> d[CCB]           -0.25 0.05 -0.35 -0.29 -0.25 -0.22 -0.15     3358     3133    1
#> d[Placebo]       -0.25 0.06 -0.36 -0.28 -0.25 -0.21 -0.14     4578     3624    1
plot(db_releff_FE, ref_line = 0)
```

![](example_diabetes_files/figure-html/diabetes_releff_FE-1.png)

``` r

(db_releff_RE <- relative_effects(db_fit_RE, trt_ref = "Diuretic"))
#>                   mean   sd  2.5%   25%   50%   75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d[Beta Blocker]  -0.07 0.09 -0.24 -0.13 -0.07 -0.01  0.10     2120     2478    1
#> d[ACE Inhibitor] -0.40 0.09 -0.58 -0.46 -0.40 -0.34 -0.24     4445     2588    1
#> d[ARB]           -0.47 0.11 -0.70 -0.54 -0.47 -0.40 -0.26     3452     2448    1
#> d[CCB]           -0.24 0.08 -0.41 -0.30 -0.24 -0.18 -0.07     4137     3318    1
#> d[Placebo]       -0.29 0.09 -0.48 -0.34 -0.29 -0.22 -0.12     4330     2852    1
plot(db_releff_RE, ref_line = 0)
```

![](example_diabetes_files/figure-html/diabetes_releff_RE-1.png)

Dias et al. ([2011](#ref-TSD2)) produce absolute predictions of the
probability of developing diabetes after three years, assuming a Normal
distribution on the baseline cloglog probability of developing diabetes
on diuretic treatment with mean -4.2 and precision 1.11. We can
replicate these results using the
[`predict()`](https://rdrr.io/r/stats/predict.html) method. We specify a
data frame of `newdata`, containing the `time` offset(s) at which to
produce predictions (here only 3 years). The `baseline` argument takes a
[`distr()`](https://dmphillippo.github.io/multinma/dev/reference/distr.md)
distribution object with which we specify the corresponding Normal
distribution on the baseline cloglog probability, and we set
`baseline_trt = "Diuretic"` to indicate that the baseline distribution
corresponds to “Diuretic” rather than the network reference “Beta
Blocker”. We set `type = "response"` to produce predicted event
probabilities (`type = "link"` would produce predicted cloglog
probabilities).

``` r

db_pred_FE <- predict(db_fit_FE, 
                      newdata = data.frame(time = 3),
                      baseline = distr(qnorm, mean = -4.2, sd = 1.11^-0.5), 
                      baseline_trt = "Diuretic",
                      type = "response")
db_pred_FE
#> ------------------------------------------------------------------ Study: New 1 ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[New 1: Beta Blocker]  0.06 0.06 0.01 0.02 0.04 0.08  0.22     4072     3918    1
#> pred[New 1: ACE Inhibitor] 0.05 0.05 0.01 0.02 0.03 0.06  0.17     4089     3919    1
#> pred[New 1: ARB]           0.04 0.04 0.00 0.02 0.03 0.05  0.16     4034     3878    1
#> pred[New 1: CCB]           0.05 0.05 0.01 0.02 0.03 0.06  0.19     4072     3919    1
#> pred[New 1: Diuretic]      0.06 0.06 0.01 0.02 0.04 0.08  0.23     4064     3878    1
#> pred[New 1: Placebo]       0.05 0.05 0.01 0.02 0.03 0.06  0.19     4065     3879    1
plot(db_pred_FE)
```

![](example_diabetes_files/figure-html/db_pred_FE-1.png)

``` r

db_pred_RE <- predict(db_fit_RE, 
                      newdata = data.frame(time = 3),
                      baseline = distr(qnorm, mean = -4.2, sd = 1.11^-0.5), 
                      baseline_trt = "Diuretic",
                      type = "response")
db_pred_RE
#> ------------------------------------------------------------------ Study: New 1 ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[New 1: Beta Blocker]  0.06 0.07 0.01 0.02 0.04 0.08  0.24     4384     3681    1
#> pred[New 1: ACE Inhibitor] 0.05 0.05 0.00 0.02 0.03 0.06  0.18     4394     3930    1
#> pred[New 1: ARB]           0.04 0.05 0.00 0.02 0.03 0.05  0.17     4365     3829    1
#> pred[New 1: CCB]           0.05 0.06 0.01 0.02 0.04 0.07  0.21     4384     3678    1
#> pred[New 1: Diuretic]      0.07 0.07 0.01 0.02 0.04 0.08  0.26     4365     3649    1
#> pred[New 1: Placebo]       0.05 0.06 0.01 0.02 0.03 0.06  0.20     4386     3780    1
plot(db_pred_RE)
```

![](example_diabetes_files/figure-html/db_pred_RE-1.png)

If the `baseline` and `newdata` arguments are omitted, predicted
probabilities will be produced for every study in the network based on
their follow-up times and estimated baseline cloglog probabilities
\mu_j:

``` r

db_pred_RE_studies <- predict(db_fit_RE, type = "response")
db_pred_RE_studies
#> ------------------------------------------------------------------- Study: AASK ---- 
#> 
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[AASK: Beta Blocker]  0.17 0.02 0.14 0.16 0.17 0.18  0.20     5901     3262    1
#> pred[AASK: ACE Inhibitor] 0.12 0.01 0.10 0.12 0.12 0.13  0.15     4607     3013    1
#> pred[AASK: ARB]           0.12 0.01 0.09 0.11 0.12 0.13  0.15     4192     2887    1
#> pred[AASK: CCB]           0.15 0.01 0.12 0.13 0.14 0.15  0.18     5481     3140    1
#> pred[AASK: Diuretic]      0.18 0.02 0.14 0.17 0.18 0.19  0.22     4467     3192    1
#> pred[AASK: Placebo]       0.14 0.02 0.11 0.13 0.14 0.15  0.17     4327     3266    1
#> 
#> ----------------------------------------------------------------- Study: ALLHAT ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[ALLHAT: Beta Blocker]  0.04 0.01 0.03 0.04 0.04 0.05  0.06     2928     2613    1
#> pred[ALLHAT: ACE Inhibitor] 0.03 0.00 0.02 0.03 0.03 0.03  0.04     4071     2634    1
#> pred[ALLHAT: ARB]           0.03 0.00 0.02 0.03 0.03 0.03  0.04     3752     2530    1
#> pred[ALLHAT: CCB]           0.04 0.00 0.03 0.03 0.04 0.04  0.05     3929     2804    1
#> pred[ALLHAT: Diuretic]      0.05 0.01 0.04 0.04 0.05 0.05  0.06     4497     2952    1
#> pred[ALLHAT: Placebo]       0.03 0.00 0.03 0.03 0.03 0.04  0.04     3797     2493    1
#> 
#> ----------------------------------------------------------------- Study: ALPINE ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[ALPINE: Beta Blocker]  0.03 0.01 0.01 0.02 0.03 0.03  0.05     6163     3107    1
#> pred[ALPINE: ACE Inhibitor] 0.02 0.01 0.01 0.01 0.02 0.02  0.03     6487     3135    1
#> pred[ALPINE: ARB]           0.02 0.01 0.01 0.01 0.02 0.02  0.03     6392     3212    1
#> pred[ALPINE: CCB]           0.02 0.01 0.01 0.02 0.02 0.03  0.04     6788     3082    1
#> pred[ALPINE: Diuretic]      0.03 0.01 0.01 0.02 0.03 0.03  0.05     6901     3139    1
#> pred[ALPINE: Placebo]       0.02 0.01 0.01 0.02 0.02 0.03  0.04     6632     3043    1
#> 
#> ----------------------------------------------------------------- Study: ANBP-2 ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[ANBP-2: Beta Blocker]  0.07 0.01 0.05 0.06 0.07 0.07  0.09     3029     2447    1
#> pred[ANBP-2: ACE Inhibitor] 0.05 0.01 0.04 0.04 0.05 0.05  0.06     4437     3017    1
#> pred[ANBP-2: ARB]           0.05 0.01 0.03 0.04 0.05 0.05  0.06     3839     2578    1
#> pred[ANBP-2: CCB]           0.06 0.01 0.04 0.05 0.06 0.06  0.07     4298     2909    1
#> pred[ANBP-2: Diuretic]      0.07 0.01 0.06 0.07 0.07 0.08  0.09     4971     3042    1
#> pred[ANBP-2: Placebo]       0.05 0.01 0.04 0.05 0.05 0.06  0.07     4431     2871    1
#> 
#> ------------------------------------------------------------------ Study: ASCOT ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[ASCOT: Beta Blocker]  0.11 0.00 0.10 0.11 0.11 0.11  0.12     6494     3085    1
#> pred[ASCOT: ACE Inhibitor] 0.08 0.01 0.07 0.08 0.08 0.09  0.10     2291     2461    1
#> pred[ASCOT: ARB]           0.08 0.01 0.06 0.07 0.08 0.08  0.09     2507     2322    1
#> pred[ASCOT: CCB]           0.10 0.01 0.08 0.09 0.10 0.10  0.11     2558     2923    1
#> pred[ASCOT: Diuretic]      0.12 0.01 0.10 0.11 0.12 0.13  0.14     2431     2940    1
#> pred[ASCOT: Placebo]       0.09 0.01 0.08 0.09 0.09 0.10  0.11     2132     2776    1
#> 
#> ------------------------------------------------------------------ Study: CAPPP ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[CAPPP: Beta Blocker]  0.07 0.00 0.07 0.07 0.07 0.08  0.08     5783     3079 1.01
#> pred[CAPPP: ACE Inhibitor] 0.05 0.00 0.05 0.05 0.05 0.06  0.06     2212     2006 1.00
#> pred[CAPPP: ARB]           0.05 0.01 0.04 0.05 0.05 0.05  0.06     2474     2354 1.00
#> pred[CAPPP: CCB]           0.06 0.00 0.05 0.06 0.06 0.07  0.07     3106     3123 1.00
#> pred[CAPPP: Diuretic]      0.08 0.01 0.07 0.08 0.08 0.08  0.10     2640     2634 1.00
#> pred[CAPPP: Placebo]       0.06 0.01 0.05 0.06 0.06 0.06  0.07     2111     2730 1.00
#> 
#> ------------------------------------------------------------------ Study: CHARM ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[CHARM: Beta Blocker]  0.09 0.01 0.07 0.08 0.09 0.10  0.12     2860     2468    1
#> pred[CHARM: ACE Inhibitor] 0.07 0.01 0.05 0.06 0.07 0.07  0.09     3556     2378    1
#> pred[CHARM: ARB]           0.06 0.01 0.05 0.06 0.06 0.07  0.08     4652     2996    1
#> pred[CHARM: CCB]           0.08 0.01 0.06 0.07 0.08 0.08  0.10     3803     2830    1
#> pred[CHARM: Diuretic]      0.10 0.02 0.07 0.09 0.10 0.10  0.13     3789     2497    1
#> pred[CHARM: Placebo]       0.07 0.01 0.06 0.07 0.07 0.08  0.09     4305     2506    1
#> 
#> ------------------------------------------------------------------ Study: DREAM ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[DREAM: Beta Blocker]  0.23 0.03 0.18 0.21 0.23 0.24  0.29     3016     2285    1
#> pred[DREAM: ACE Inhibitor] 0.17 0.02 0.13 0.16 0.17 0.18  0.21     4980     2874    1
#> pred[DREAM: ARB]           0.16 0.02 0.12 0.15 0.16 0.17  0.21     4220     2383    1
#> pred[DREAM: CCB]           0.20 0.03 0.15 0.18 0.19 0.21  0.25     4383     2866    1
#> pred[DREAM: Diuretic]      0.24 0.03 0.19 0.22 0.24 0.26  0.31     4958     2492    1
#> pred[DREAM: Placebo]       0.19 0.02 0.15 0.17 0.19 0.20  0.23     5096     2976    1
#> 
#> ------------------------------------------------------------------- Study: EWPH ---- 
#> 
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[EWPH: Beta Blocker]  0.06 0.01 0.04 0.05 0.06 0.07  0.09     3960     2586    1
#> pred[EWPH: ACE Inhibitor] 0.05 0.01 0.03 0.04 0.04 0.05  0.06     5446     3085    1
#> pred[EWPH: ARB]           0.04 0.01 0.03 0.04 0.04 0.05  0.06     5099     2476    1
#> pred[EWPH: CCB]           0.05 0.01 0.04 0.05 0.05 0.06  0.08     5268     2576    1
#> pred[EWPH: Diuretic]      0.07 0.01 0.05 0.06 0.07 0.07  0.09     5430     3076    1
#> pred[EWPH: Placebo]       0.05 0.01 0.03 0.04 0.05 0.06  0.07     5815     2613    1
#> 
#> ------------------------------------------------------------------ Study: FEVER ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[FEVER: Beta Blocker]  0.04 0.01 0.03 0.04 0.04 0.04  0.05     2955     2712    1
#> pred[FEVER: ACE Inhibitor] 0.03 0.00 0.02 0.03 0.03 0.03  0.04     4989     2866    1
#> pred[FEVER: ARB]           0.03 0.00 0.02 0.03 0.03 0.03  0.04     4534     2710    1
#> pred[FEVER: CCB]           0.04 0.00 0.03 0.03 0.03 0.04  0.05     4676     3035    1
#> pred[FEVER: Diuretic]      0.04 0.01 0.03 0.04 0.04 0.05  0.06     4621     2864    1
#> pred[FEVER: Placebo]       0.03 0.00 0.03 0.03 0.03 0.04  0.04     5193     3183    1
#> 
#> ----------------------------------------------------------------- Study: HAPPHY ---- 
#> 
#>                             mean sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[HAPPHY: Beta Blocker]  0.02  0 0.02 0.02 0.02 0.03  0.03     6053     2917    1
#> pred[HAPPHY: ACE Inhibitor] 0.02  0 0.01 0.02 0.02 0.02  0.02     4802     3209    1
#> pred[HAPPHY: ARB]           0.02  0 0.01 0.02 0.02 0.02  0.02     4102     2806    1
#> pred[HAPPHY: CCB]           0.02  0 0.02 0.02 0.02 0.02  0.03     5170     3236    1
#> pred[HAPPHY: Diuretic]      0.03  0 0.02 0.02 0.03 0.03  0.03     4284     3215    1
#> pred[HAPPHY: Placebo]       0.02  0 0.02 0.02 0.02 0.02  0.03     4507     2599    1
#> 
#> ------------------------------------------------------------------- Study: HOPE ---- 
#> 
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[HOPE: Beta Blocker]  0.06 0.01 0.04 0.05 0.06 0.06  0.08     3180     2550    1
#> pred[HOPE: ACE Inhibitor] 0.04 0.01 0.03 0.04 0.04 0.05  0.06     4869     2828    1
#> pred[HOPE: ARB]           0.04 0.01 0.03 0.04 0.04 0.04  0.05     4281     2769    1
#> pred[HOPE: CCB]           0.05 0.01 0.04 0.04 0.05 0.05  0.07     4548     2830    1
#> pred[HOPE: Diuretic]      0.06 0.01 0.05 0.06 0.06 0.07  0.08     4582     2978    1
#> pred[HOPE: Placebo]       0.05 0.01 0.04 0.04 0.05 0.05  0.06     5631     3059    1
#> 
#> ---------------------------------------------------------------- Study: INSIGHT ---- 
#> 
#>                              mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[INSIGHT: Beta Blocker]  0.07 0.01 0.05 0.06 0.06 0.07  0.09     2860     2638    1
#> pred[INSIGHT: ACE Inhibitor] 0.05 0.01 0.03 0.04 0.05 0.05  0.06     3756     2716    1
#> pred[INSIGHT: ARB]           0.04 0.01 0.03 0.04 0.04 0.05  0.06     3836     2376    1
#> pred[INSIGHT: CCB]           0.06 0.01 0.04 0.05 0.05 0.06  0.07     4002     2676    1
#> pred[INSIGHT: Diuretic]      0.07 0.01 0.05 0.06 0.07 0.08  0.09     4344     2947    1
#> pred[INSIGHT: Placebo]       0.05 0.01 0.04 0.05 0.05 0.06  0.07     4062     2693    1
#> 
#> ----------------------------------------------------------------- Study: INVEST ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[INVEST: Beta Blocker]  0.08 0.00 0.08 0.08 0.08 0.08  0.09     7301     2876    1
#> pred[INVEST: ACE Inhibitor] 0.06 0.00 0.05 0.06 0.06 0.06  0.07     2417     2558    1
#> pred[INVEST: ARB]           0.06 0.01 0.05 0.05 0.06 0.06  0.07     2495     2453    1
#> pred[INVEST: CCB]           0.07 0.00 0.06 0.07 0.07 0.07  0.08     2698     2785    1
#> pred[INVEST: Diuretic]      0.09 0.01 0.07 0.08 0.09 0.09  0.11     2458     2573    1
#> pred[INVEST: Placebo]       0.07 0.01 0.06 0.06 0.07 0.07  0.08     2177     2812    1
#> 
#> ------------------------------------------------------------------- Study: LIFE ---- 
#> 
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[LIFE: Beta Blocker]  0.08 0.00 0.07 0.08 0.08 0.08  0.09     7399     3174    1
#> pred[LIFE: ACE Inhibitor] 0.06 0.01 0.05 0.06 0.06 0.06  0.07     2681     2387    1
#> pred[LIFE: ARB]           0.06 0.01 0.05 0.05 0.06 0.06  0.07     2532     2242    1
#> pred[LIFE: CCB]           0.07 0.01 0.06 0.07 0.07 0.07  0.08     3479     2916    1
#> pred[LIFE: Diuretic]      0.09 0.01 0.07 0.08 0.09 0.09  0.10     2901     2432    1
#> pred[LIFE: Placebo]       0.07 0.01 0.05 0.06 0.07 0.07  0.08     2411     2599    1
#> 
#> ------------------------------------------------------------------ Study: MRC-E ---- 
#> 
#>                            mean sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[MRC-E: Beta Blocker]  0.03  0 0.02 0.03 0.03 0.03  0.04     4193     3016    1
#> pred[MRC-E: ACE Inhibitor] 0.02  0 0.02 0.02 0.02 0.02  0.03     5088     3141    1
#> pred[MRC-E: ARB]           0.02  0 0.01 0.02 0.02 0.02  0.03     5211     2894    1
#> pred[MRC-E: CCB]           0.03  0 0.02 0.02 0.02 0.03  0.03     5146     3267    1
#> pred[MRC-E: Diuretic]      0.03  0 0.02 0.03 0.03 0.03  0.04     4272     3362    1
#> pred[MRC-E: Placebo]       0.02  0 0.02 0.02 0.02 0.03  0.03     4904     3238    1
#> 
#> ----------------------------------------------------------------- Study: NORDIL ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[NORDIL: Beta Blocker]  0.05 0.00 0.04 0.05 0.05 0.05  0.06     6985     3046    1
#> pred[NORDIL: ACE Inhibitor] 0.04 0.00 0.03 0.03 0.04 0.04  0.04     2902     2488    1
#> pred[NORDIL: ARB]           0.03 0.00 0.03 0.03 0.03 0.04  0.04     2926     2654    1
#> pred[NORDIL: CCB]           0.04 0.00 0.04 0.04 0.04 0.04  0.05     3285     2874    1
#> pred[NORDIL: Diuretic]      0.05 0.01 0.04 0.05 0.05 0.06  0.06     2823     2667    1
#> pred[NORDIL: Placebo]       0.04 0.00 0.03 0.04 0.04 0.04  0.05     2613     2971    1
#> 
#> ------------------------------------------------------------------ Study: PEACE ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[PEACE: Beta Blocker]  0.14 0.02 0.11 0.13 0.14 0.15  0.18     2855     2335    1
#> pred[PEACE: ACE Inhibitor] 0.10 0.01 0.08 0.09 0.10 0.11  0.13     4699     2481    1
#> pred[PEACE: ARB]           0.09 0.01 0.07 0.09 0.09 0.10  0.13     4311     2738    1
#> pred[PEACE: CCB]           0.12 0.02 0.09 0.11 0.12 0.13  0.15     4017     2231    1
#> pred[PEACE: Diuretic]      0.15 0.02 0.11 0.13 0.15 0.16  0.19     3954     2423    1
#> pred[PEACE: Placebo]       0.11 0.01 0.09 0.10 0.11 0.12  0.14     5177     2685    1
#> 
#> ------------------------------------------------------------------ Study: SCOPE ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[SCOPE: Beta Blocker]  0.06 0.01 0.05 0.06 0.06 0.07  0.09     3151     2729    1
#> pred[SCOPE: ACE Inhibitor] 0.05 0.01 0.03 0.04 0.05 0.05  0.06     4625     3115    1
#> pred[SCOPE: ARB]           0.04 0.01 0.03 0.04 0.04 0.05  0.06     4936     3072    1
#> pred[SCOPE: CCB]           0.05 0.01 0.04 0.05 0.05 0.06  0.07     4254     2781    1
#> pred[SCOPE: Diuretic]      0.07 0.01 0.05 0.06 0.07 0.08  0.09     4086     3162    1
#> pred[SCOPE: Placebo]       0.05 0.01 0.04 0.05 0.05 0.06  0.07     4957     2783    1
#> 
#> ------------------------------------------------------------------- Study: SHEP ---- 
#> 
#>                           mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[SHEP: Beta Blocker]  0.09 0.01 0.06 0.08 0.09 0.09  0.11     2980     2717    1
#> pred[SHEP: ACE Inhibitor] 0.06 0.01 0.05 0.06 0.06 0.07  0.08     4707     3090    1
#> pred[SHEP: ARB]           0.06 0.01 0.04 0.05 0.06 0.06  0.08     4486     3032    1
#> pred[SHEP: CCB]           0.07 0.01 0.05 0.07 0.07 0.08  0.10     4647     3027    1
#> pred[SHEP: Diuretic]      0.09 0.01 0.07 0.08 0.09 0.10  0.12     5141     3014    1
#> pred[SHEP: Placebo]       0.07 0.01 0.05 0.06 0.07 0.08  0.09     4734     2849    1
#> 
#> ----------------------------------------------------------------- Study: STOP-2 ---- 
#> 
#>                             mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[STOP-2: Beta Blocker]  0.05 0.00 0.05 0.05 0.05 0.06  0.06     4797     2959    1
#> pred[STOP-2: ACE Inhibitor] 0.04 0.00 0.03 0.04 0.04 0.04  0.05     2906     2629    1
#> pred[STOP-2: ARB]           0.04 0.00 0.03 0.03 0.04 0.04  0.04     2810     2496    1
#> pred[STOP-2: CCB]           0.05 0.00 0.04 0.04 0.05 0.05  0.05     3670     2662    1
#> pred[STOP-2: Diuretic]      0.06 0.01 0.05 0.05 0.06 0.06  0.07     3365     2682    1
#> pred[STOP-2: Placebo]       0.04 0.00 0.03 0.04 0.04 0.05  0.05     2711     2868    1
#> 
#> ------------------------------------------------------------------ Study: VALUE ---- 
#> 
#>                            mean   sd 2.5%  25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> pred[VALUE: Beta Blocker]  0.20 0.02 0.15 0.18 0.19 0.21  0.25     2933     2458    1
#> pred[VALUE: ACE Inhibitor] 0.15 0.02 0.11 0.13 0.14 0.16  0.19     3770     2582    1
#> pred[VALUE: ARB]           0.14 0.02 0.10 0.13 0.14 0.15  0.17     4152     2795    1
#> pred[VALUE: CCB]           0.17 0.02 0.13 0.16 0.17 0.18  0.21     4056     2601    1
#> pred[VALUE: Diuretic]      0.21 0.03 0.16 0.19 0.21 0.22  0.27     3866     2659    1
#> pred[VALUE: Placebo]       0.16 0.02 0.12 0.15 0.16 0.17  0.21     3909     2445    1
plot(db_pred_RE_studies)
```

![](example_diabetes_files/figure-html/db_pred_RE_all-1.png)

We can also produce treatment rankings, rank probabilities, and
cumulative rank probabilities.

``` r

(db_ranks <- posterior_ranks(db_fit_RE))
#>                     mean   sd 2.5% 25% 50% 75% 97.5% Bulk_ESS Tail_ESS Rhat
#> rank[Beta Blocker]  5.20 0.42    5   5   5   5     6     2226       NA    1
#> rank[ACE Inhibitor] 1.84 0.53    1   2   2   2     3     3407     3760    1
#> rank[ARB]           1.27 0.53    1   1   1   1     3     3321     3017    1
#> rank[CCB]           3.70 0.52    3   3   4   4     4     3958     3154    1
#> rank[Diuretic]      5.79 0.42    5   6   6   6     6     2438       NA    1
#> rank[Placebo]       3.20 0.59    2   3   3   4     4     4049     3313    1
plot(db_ranks)
```

![](example_diabetes_files/figure-html/diabetes_ranks-1.png)

``` r

(db_rankprobs <- posterior_rank_probs(db_fit_RE))
#>                  p_rank[1] p_rank[2] p_rank[3] p_rank[4] p_rank[5] p_rank[6]
#> d[Beta Blocker]       0.00      0.00      0.00      0.01      0.79       0.2
#> d[ACE Inhibitor]      0.22      0.71      0.06      0.00      0.00       0.0
#> d[ARB]                0.77      0.20      0.03      0.00      0.00       0.0
#> d[CCB]                0.00      0.02      0.27      0.71      0.01       0.0
#> d[Diuretic]           0.00      0.00      0.00      0.00      0.20       0.8
#> d[Placebo]            0.01      0.07      0.65      0.27      0.00       0.0
plot(db_rankprobs)
```

![](example_diabetes_files/figure-html/diabetes_rankprobs-1.png)

``` r

(db_cumrankprobs <- posterior_rank_probs(db_fit_RE, cumulative = TRUE))
#>                  p_rank[1] p_rank[2] p_rank[3] p_rank[4] p_rank[5] p_rank[6]
#> d[Beta Blocker]       0.00      0.00      0.00      0.01       0.8         1
#> d[ACE Inhibitor]      0.22      0.94      1.00      1.00       1.0         1
#> d[ARB]                0.77      0.97      1.00      1.00       1.0         1
#> d[CCB]                0.00      0.02      0.28      0.99       1.0         1
#> d[Diuretic]           0.00      0.00      0.00      0.00       0.2         1
#> d[Placebo]            0.01      0.08      0.72      1.00       1.0         1
plot(db_cumrankprobs)
```

![](example_diabetes_files/figure-html/diabetes_cumrankprobs-1.png)

## References

Dias, S., N. J. Welton, A. J. Sutton, and A. E. Ades. 2011. *NICE DSU
Technical Support Document 2: A Generalised Linear Modelling Framework
for Pair-Wise and Network Meta-Analysis of Randomised Controlled
Trials*. National Institute for Health and Care Excellence.
<https://sheffield.ac.uk/nice-dsu>.

Elliott, W. J., and P. M. Meyer. 2007. “Incident Diabetes in Clinical
Trials of Antihypertensive Drugs: A Network Meta-Analysis.” *The Lancet*
369 (9557): 201–7. <https://doi.org/10.1016/s0140-6736(07)60108-1>.

[^1]: The gain in efficiency here from using “Beta Blocker” as the
    network reference treatment instead of “Diuretic” is considerable -
    around 4-8 times, in terms of effective samples per second. The
    functions in this package will always attempt to choose a default
    network reference treatment that maximises computational efficiency
    and stability. If you have chosen an alternative network reference
    treatment and the model runs very slowly or has low effective sample
    size, this is a likely cause.
