# Network meta-analysis models

The `nma` function fits network meta-analysis and (multilevel) network
meta-regression models in Stan.

## Usage

``` r
nma(
  network,
  consistency = c("consistency", "ume", "nodesplit"),
  trt_effects = c("fixed", "random"),
  regression = NULL,
  class_interactions = c("common", "exchangeable", "independent"),
  class_effects = c("independent", "common", "exchangeable"),
  class_sd = c("independent", "common"),
  likelihood = NULL,
  link = NULL,
  ...,
  nodesplit = get_nodesplits(network, include_consistency = TRUE),
  prior_intercept = .default(normal(scale = 100)),
  prior_trt = .default(normal(scale = 10)),
  prior_het = .default(half_normal(scale = 5)),
  prior_het_type = c("sd", "var", "prec"),
  prior_reg = .default(normal(scale = 10)),
  prior_aux = .default(),
  prior_aux_reg = .default(),
  prior_class_mean = .default(normal(scale = 10)),
  prior_class_sd = .default(half_normal(scale = 5)),
  aux_by = NULL,
  aux_regression = NULL,
  QR = FALSE,
  center = TRUE,
  adapt_delta = NULL,
  int_thin = 0,
  int_check = TRUE,
  mspline_degree = 3,
  n_knots = 7,
  knots = NULL,
  mspline_basis = NULL
)
```

## Arguments

- network:

  An `nma_data` object, as created by the functions `set_*()`,
  [`combine_network()`](https://dmphillippo.github.io/multinma/reference/combine_network.md),
  or
  [`add_integration()`](https://dmphillippo.github.io/multinma/reference/add_integration.md)

- consistency:

  Character string specifying the type of (in)consistency model to fit,
  either `"consistency"`, `"ume"`, or `"nodesplit"`

- trt_effects:

  Character string specifying either `"fixed"` or `"random"` effects

- regression:

  A one-sided model formula, specifying the prognostic and
  effect-modifying terms for a regression model. Any references to
  treatment should use the `.trt` special variable, for example
  specifying effect modifier interactions as `variable:.trt` (see
  details). The special variable `.mu` can be used for baseline risk
  meta-regression: `~.mu:.trt`.

- class_interactions:

  Character string specifying whether effect modifier interactions are
  specified as `"common"`, `"exchangeable"`, or `"independent"`.

- class_effects:

  Character string specifying a model for treatment class effects,
  either `"independent"` (the default), `"exchangeable"`, or `"common"`.

- class_sd:

  Character string specifying whether the class standard deviations in a
  class effects model should be `"independent"` (i.e. separate for each
  class, the default), or `"common"` (i.e. shared across all classes).
  Alternatively this can be a list of character vectors, each of which
  describe a set classes for which to share a common class SD; any list
  names will be used to name the output parameters, otherwise the name
  will be taken from the first class in each set.

- likelihood:

  Character string specifying a likelihood, if unspecified will be
  inferred from the data (see details)

- link:

  Character string specifying a link function, if unspecified will
  default to the canonical link (see details)

- ...:

  Further arguments passed to `sampling()`, such as `iter`, `chains`,
  `cores`, etc.

- nodesplit:

  For `consistency = "nodesplit"`, the comparison(s) to split in the
  node-splitting model(s). Either a length 2 vector giving the
  treatments in a single comparison, or a 2 column data frame listing
  multiple treatment comparisons to split in turn. By default, all
  possible comparisons will be chosen (see
  [`get_nodesplits()`](https://dmphillippo.github.io/multinma/reference/get_nodesplits.md)).

- prior_intercept:

  Specification of prior distribution for the intercept

- prior_trt:

  Specification of prior distribution for the treatment effects

- prior_het:

  Specification of prior distribution for the heterogeneity (if
  `trt_effects = "random"`)

- prior_het_type:

  Character string specifying whether the prior distribution `prior_het`
  is placed on the heterogeneity standard deviation \\\tau\\ (`"sd"`,
  the default), variance \\\tau^2\\ (`"var"`), or precision \\1/\tau^2\\
  (`"prec"`).

- prior_reg:

  Specification of prior distribution for the regression coefficients
  (if `regression` formula specified)

- prior_aux:

  Specification of prior distribution for the auxiliary parameter, if
  applicable (see details). For `likelihood = "gengamma"` this should be
  a list of prior distributions with elements `sigma` and `k`.

- prior_aux_reg:

  Specification of prior distribution for the auxiliary regression
  parameters, if `aux_regression` is specified (see details).

- prior_class_mean:

  Specification of prior distribution for the treatment class means (if
  `class_effects = "exchangeable"`).

- prior_class_sd:

  Specification of prior distribution for the treatment class standard
  deviations (if `class_effects = "exchangeable"`).

- aux_by:

  Vector of variable names listing the variables to stratify the
  auxiliary parameters by. Currently only used for survival models, see
  details. Cannot be used with `aux_regression`.

- aux_regression:

  A one-sided model formula giving a regression model for the auxiliary
  parameters. Currently only used for survival models, see details.
  Cannot be used with `aux_by`.

- QR:

  Logical scalar (default `FALSE`), whether to apply a QR decomposition
  to the model design matrix

- center:

  Logical scalar (default `TRUE`), whether to center the (numeric)
  regression terms about the overall means

- adapt_delta:

  See
  [adapt_delta](https://dmphillippo.github.io/multinma/reference/adapt_delta.md)
  for details

- int_thin:

  A single integer value, the thinning factor for returning cumulative
  estimates of integration error. Saving cumulative estimates is
  disabled by `int_thin = 0`, which is the default.

- int_check:

  Logical, check sufficient accuracy of numerical integration by fitting
  half of the chains with `n_int/2`? When `TRUE`, `Rhat` and `n_eff`
  diagnostic warnings will be given if numerical integration has not
  sufficiently converged (suggesting increasing `n_int` in
  [`add_integration()`](https://dmphillippo.github.io/multinma/reference/add_integration.md)).
  Default `TRUE`, but disabled (`FALSE`) when `int_thin > 0`.

- mspline_degree:

  Non-negative integer giving the degree of the M-spline polynomial for
  `likelihood = "mspline"`. Piecewise exponential hazards
  (`likelihood = "pexp"`) are a special case with `mspline_degree = 0`.

- n_knots:

  For `mspline` and `pexp` likelihoods, a non-negative integer giving
  the number of internal knots for partitioning the baseline hazard into
  intervals. The knot locations within each study will be determined by
  the corresponding quantiles of the observed event times, plus boundary
  knots at the earliest entry time (0 with no delayed entry) and the
  maximum event/censoring time. For example, with `n_knots = 3`, the
  internal knot locations will be at the 25%, 50%, and 75% quantiles of
  the observed event times. The default is `n_knots = 7`; overfitting is
  avoided by shrinking towards a constant hazard with a random walk
  prior (see details). If `aux_regression` is specified then a single
  set of knot locations will be calculated across all studies in the
  network. See
  [`make_knots()`](https://dmphillippo.github.io/multinma/reference/make_knots.md)
  for more details on the knot positioning algorithms. Ignored when
  `knots` is specified.

- knots:

  For `mspline` and `pexp` likelihoods, a named list of numeric vectors
  of knot locations (including boundary knots) for each of the studies
  in the network. Currently, each vector must have the same length (i.e.
  each study must use the same number of knots). Alternatively, a single
  numeric vector of knot locations can be provided which will be shared
  across all studies in the network. If unspecified (the default), the
  knots will be chosen based on `n_knots` as described above. If
  `aux_regression` is specified then `knots` should be a single numeric
  vector of knot locations which will be shared across all studies in
  the network.
  [`make_knots()`](https://dmphillippo.github.io/multinma/reference/make_knots.md)
  can be used to help specify `knots` directly, or to investigate knot
  placement prior to model fitting.

- mspline_basis:

  Instead of specifying `mspline_degree` and `n_knots` or `knots`, a
  named list of M-spline bases (one for each study) can be provided with
  `mspline_basis` which will be used directly. In this case, all other
  M-spline options will be ignored.

## Value

`nma()` returns a
[stan_nma](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
object, except when `consistency = "nodesplit"` when a
[nma_nodesplit](https://dmphillippo.github.io/multinma/reference/nma_nodesplit-class.md)
or
[nma_nodesplit_df](https://dmphillippo.github.io/multinma/reference/nma_nodesplit-class.md)
object is returned. `nma.fit()` returns a
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
object.

## Details

When specifying a model formula in the `regression` argument, the usual
formula syntax is available (as interpreted by
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)). The only
additional requirement here is that the special variable `.trt` should
be used to refer to treatment. For example, effect modifier interactions
should be specified as `variable:.trt`. Prognostic (main) effects and
interactions can be included together compactly as `variable*.trt`,
which expands to `variable + variable:.trt` (plus `.trt`, which is
already in the NMA model).

For the advanced user, the additional specials `.study` and `.trtclass`
are also available, and refer to studies and (if specified) treatment
classes respectively. When node-splitting models are fitted
(`consistency = "nodesplit"`) the special `.omega` is available,
indicating the arms to which the node-splitting inconsistency factor is
added.

See
[`?priors`](https://dmphillippo.github.io/multinma/reference/priors.md)
for details on prior specification. Default prior distributions are
available, but may not be appropriate for the particular setting and
will raise a warning if used. No attempt is made to tailor these
defaults to the data provided. Please consider appropriate prior
distributions for the particular setting, accounting for the scales of
outcomes and covariates, etc. The function
[`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/reference/plot_prior_posterior.md)
may be useful in examining the influence of the chosen prior
distributions on the posterior distributions, and the
[`summary()`](https://dmphillippo.github.io/multinma/reference/summary.nma_prior.md)
method for `nma_prior` objects prints prior intervals.

## Likelihoods and link functions

Currently, the following likelihoods and link functions are supported
for each data type:

|                |                                                                                                                                            |                              |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------|------------------------------|
| **Data type**  | **Likelihood**                                                                                                                             | **Link function**            |
| **Binary**     | `bernoulli`, `bernoulli2`                                                                                                                  | `logit`, `probit`, `cloglog` |
| **Count**      | `binomial`, `binomial2`                                                                                                                    | `logit`, `probit`, `cloglog` |
| **Rate**       | `poisson`                                                                                                                                  | `log`                        |
| **Continuous** | `normal`                                                                                                                                   | `identity`, `log`            |
| **Ordered**    | `ordered`                                                                                                                                  | `logit`, `probit`, `cloglog` |
| **Survival**   | `exponential`, `weibull`, `gompertz`, `exponential-aft`, `weibull-aft`, `lognormal`, `loglogistic`, `gamma`, `gengamma`, `mspline`, `pexp` | `log`                        |

The `bernoulli2` and `binomial2` likelihoods correspond to a
two-parameter Binomial likelihood for arm-based AgD, which more closely
matches the underlying Poisson Binomial distribution for the summarised
aggregate outcomes in a ML-NMR model than the typical (one parameter)
Binomial distribution (see Phillippo et al. 2020) .

When a `cloglog` link is used, including an offset for log follow-up
time (i.e. `regression = ~offset(log(time))`) results in a model on the
log hazard (see Dias et al. 2011) .

For survival data, all accelerated failure time models
(`exponential-aft`, `weibull-aft`, `lognormal`, `loglogistic`, `gamma`,
`gengamma`) are parameterised so that the treatment effects and any
regression parameters are log Survival Time Ratios (i.e. a coefficient
of \\\log(2)\\ means that the treatment or covariate is associated with
a doubling of expected survival time). These can be converted to log
Acceleration Factors using the relation \\\log(\mathrm{AF}) =
-\log(\mathrm{STR})\\ (or equivalently \\\mathrm{AF} =
1/\mathrm{STR}\\).

Further details on each likelihood and link function are given by Dias
et al. (2011) . Details of survival models are given by Phillippo et al.
(2025) .

## Auxiliary parameters

Auxiliary parameters are only present in the following models.

### Normal likelihood with IPD

When a Normal likelihood is fitted to IPD, the auxiliary parameters are
the arm-level standard deviations \\\sigma\_{jk}\\ on treatment \\k\\ in
study \\j\\.

### Ordered multinomial likelihood

When fitting a model to \\M\\ ordered outcomes, the auxiliary parameters
are the latent cutoffs between each category, \\c_0 \< c_1 \< \dots \<
c_M\\. Only \\c_2\\ to \\c\_{M-1}\\ are estimated; we fix \\c_0 =
-\infty\\, \\c_1 = 0\\, and \\c_M = \infty\\. When specifying priors for
these latent cutoffs, we choose to specify priors on the *differences*
\\c\_{m+1} - c_m\\. Stan automatically truncates any priors so that the
ordering constraints are satisfied.

### Survival (time-to-event) likelihoods

All survival likelihoods except the `exponential` and `exponential-aft`
likelihoods have auxiliary parameters. These are typically
study-specific shape parameters \\\gamma_j\>0\\, except for the
`lognormal` likelihood where the auxiliary parameters are study-specific
standard deviations on the log scale \\\sigma_j\>0\\.

The `gengamma` likelihood has two sets of auxiliary parameters,
study-specific scale parameters \\\sigma_j\>0\\ and shape parameters
\\k_j\\, following the parameterisation of Lawless (1980) , which
permits a range of behaviours for the baseline hazard including
increasing, decreasing, bathtub and arc-shaped hazards. This
parameterisation is related to that discussed by Cox et al. (2007) and
implemented in the `flexsurv` package with \\Q = k^{-0.5}\\. The
parameterisation used here effectively bounds the shape parameter \\k\\
away from numerical instabilities as \\k \rightarrow \infty\\ (i.e. away
from \\Q \rightarrow 0\\, the log-Normal distribution) via the prior
distribution. Implicitly, this parameterisation is restricted to \\Q \>
0\\ and so certain survival distributions like the inverse-Gamma and
inverse-Weibull are not part of the parameter space; however, \\Q \> 0\\
still encompasses the other survival distributions implemented in this
package.

For the `mspline` and `pexp` likelihoods, the auxiliary parameters are
the spline coefficients for each study. These form a unit simplex (i.e.
lie between 0 and 1, and sum to 1), and are given a random walk prior
distribution (Phillippo et al. 2025) . `prior_aux` specifies the
hyperprior on the random walk standard deviation \\\sigma\\ which
controls the level of smoothing of the baseline hazard, with \\\sigma =
0\\ corresponding to a constant baseline hazard.

The auxiliary parameters can be stratified by additional factors through
the `aux_by` argument. For example, to allow the shape of the baseline
hazard to vary between treatment arms as well as studies, use
`aux_by = c(".study", ".trt")`. (Technically, `.study` is always
included in the stratification even if omitted from `aux_by`, but we
choose here to make the stratification explicit.) This is a common way
of relaxing the proportional hazards assumption. The default is
equivalent to `aux_by = ".study"` which stratifies the auxiliary
parameters by study, as described above.

A regression model may be specified on the auxiliary parameters using
`aux_regression`. This is useful if we wish to model departures from
non-proportionality, rather than allowing the baseline hazards to be
completely independent using `aux_by`. This is necessary if absolute
predictions (e.g. survival curves) are required in a population for
unobserved combinations of covariates; for example, if `aux_by = .trt`
then absolute predictions may only be produced for the observed
treatment arms in each study population, whereas if
`aux_regression = ~.trt` then absolute predictions can be produced for
all treatments in any population. For `mspline` and `pexp` likelihoods,
the regression coefficients are smoothed over time using a random walk
prior to avoid overfitting (Phillippo et al. 2025) : `prior_aux_reg`
specifies the hyperprior for the random walk standard deviation. For
other parametric likelihoods, `prior_aux_reg` specifies the prior for
the auxiliary regression coefficients.

## References

Cox C, Chu H, Schneider MF, Muñoz A (2007). “Parametric survival
analysis and taxonomy of hazard functions for the generalized gamma
distribution.” *Statistics in Medicine*, **26**(23), 4352–4374.
[doi:10.1002/sim.2836](https://doi.org/10.1002/sim.2836) .  
  
Dias S, Welton NJ, Sutton AJ, Ades AE (2011). “NICE DSU Technical
Support Document 2: A generalised linear modelling framework for
pair-wise and network meta-analysis of randomised controlled trials.”
National Institute for Health and Care Excellence.
<https://sheffield.ac.uk/nice-dsu>.  
  
Lawless JF (1980). “Inference in the Generalized Gamma and Log Gamma
Distributions.” *Technometrics*, **22**(3), 409–419.
[doi:10.1080/00401706.1980.10486173](https://doi.org/10.1080/00401706.1980.10486173)
.  
  
Phillippo DM, Dias S, Ades AE, Belger M, Brnabic A, Schacht A, Saure D,
Kadziola Z, Welton NJ (2020). “Multilevel Network Meta-Regression for
population-adjusted treatment comparisons.” *Journal of the Royal
Statistical Society: Series A (Statistics in Society)*, **183**(3),
1189–1210. [doi:10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579)
.  
  
Phillippo DM, Dias S, Ades AE, Welton NJ (2025). “Multilevel network
meta-regression for general likelihoods: synthesis of individual and
aggregate data with applications to survival analysis.” *Journal of the
Royal Statistical Society Series A: Statistics in Society*.
[doi:10.1093/jrsssa/qnaf169](https://doi.org/10.1093/jrsssa/qnaf169) .  
  
Phillippo DM, Sadek A, Pedder H, Welton NJ (2025). “Network
Meta-Analysis of survival outcomes with non-proportional hazards using
flexible M-splines.” *arXiv*.
[doi:10.48550/ARXIV.2509.10383](https://doi.org/10.48550/ARXIV.2509.10383)
, 2509.10383.

## Examples

``` r
## Smoking cessation NMA
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
#> d[Group counselling]          0.84    0.00 0.18     0.51     0.72     0.84
#> d[Individual counselling]     0.76    0.00 0.06     0.65     0.73     0.76
#> d[Self-help]                  0.22    0.00 0.13    -0.04     0.14     0.22
#> lp__                      -5859.19    0.09 3.69 -5867.22 -5861.48 -5858.86
#>                                75%    97.5% n_eff Rhat
#> d[Group counselling]          0.96     1.19  2025    1
#> d[Individual counselling]     0.80     0.88  1998    1
#> d[Self-help]                  0.31     0.46  2696    1
#> lp__                      -5856.58 -5852.87  1699    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 17 08:15:38 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }

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
#> d[Group counselling]          1.10    0.01 0.44     0.27     0.80     1.09
#> d[Individual counselling]     0.84    0.01 0.25     0.38     0.68     0.84
#> d[Self-help]                  0.49    0.01 0.40    -0.29     0.23     0.48
#> lp__                      -5768.32    0.20 6.59 -5782.07 -5772.64 -5768.00
#> tau                           0.83    0.01 0.18     0.55     0.71     0.81
#>                                75%    97.5% n_eff Rhat
#> d[Group counselling]          1.39     1.96  1833    1
#> d[Individual counselling]     1.00     1.35  1170    1
#> d[Self-help]                  0.75     1.31  1654    1
#> lp__                      -5763.68 -5756.83  1081    1
#> tau                           0.93     1.26  1169    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 17 08:15:44 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }

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
#> d[Group counselling vs. No intervention]            1.13    0.02 0.80    -0.33
#> d[Individual counselling vs. No intervention]       0.91    0.01 0.28     0.40
#> d[Self-help vs. No intervention]                    0.33    0.01 0.58    -0.80
#> d[Individual counselling vs. Group counselling]    -0.31    0.01 0.61    -1.52
#> d[Self-help vs. Group counselling]                 -0.61    0.01 0.70    -2.03
#> d[Self-help vs. Individual counselling]             0.13    0.02 1.03    -1.93
#> lp__                                            -5765.69    0.21 6.37 -5778.98
#> tau                                                 0.92    0.01 0.22     0.58
#>                                                      25%      50%      75%
#> d[Group counselling vs. No intervention]            0.60     1.10     1.62
#> d[Individual counselling vs. No intervention]       0.72     0.89     1.08
#> d[Self-help vs. No intervention]                   -0.05     0.33     0.70
#> d[Individual counselling vs. Group counselling]    -0.70    -0.32     0.08
#> d[Self-help vs. Group counselling]                 -1.07    -0.62    -0.15
#> d[Self-help vs. Individual counselling]            -0.55     0.14     0.81
#> lp__                                            -5769.77 -5765.33 -5761.22
#> tau                                                 0.77     0.90     1.04
#>                                                    97.5% n_eff Rhat
#> d[Group counselling vs. No intervention]            2.85  2455    1
#> d[Individual counselling vs. No intervention]       1.49  1249    1
#> d[Self-help vs. No intervention]                    1.45  2489    1
#> d[Individual counselling vs. Group counselling]     0.88  2588    1
#> d[Self-help vs. Group counselling]                  0.73  2915    1
#> d[Self-help vs. Individual counselling]             2.12  3099    1
#> lp__                                            -5754.24   929    1
#> tau                                                 1.43  1050    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 17 08:15:51 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }

# \donttest{
# Fitting all possible node-splitting models
smk_fit_RE_nodesplit <- nma(smk_net, 
                            consistency = "nodesplit",
                            trt_effects = "random",
                            prior_intercept = normal(scale = 100),
                            prior_trt = normal(scale = 100),
                            prior_het = normal(scale = 5))
#> Fitting model 1 of 7, node-split: Group counselling vs. No intervention
#> Fitting model 2 of 7, node-split: Individual counselling vs. No intervention
#> Fitting model 3 of 7, node-split: Self-help vs. No intervention
#> Fitting model 4 of 7, node-split: Individual counselling vs. Group counselling
#> Fitting model 5 of 7, node-split: Self-help vs. Group counselling
#> Fitting model 6 of 7, node-split: Self-help vs. Individual counselling
#> Fitting model 7 of 7, consistency model
# }

# \donttest{
# Summarise the node-splitting results
summary(smk_fit_RE_nodesplit)
#> Node-splitting models fitted for 6 comparisons.
#> 
#> ------------------------------ Node-split Group counselling vs. No intervention ---- 
#> 
#>                  mean   sd  2.5%   25%   50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net            1.12 0.44  0.29  0.83  1.11 1.41  2.00     2316     2639    1
#> d_dir            1.07 0.76 -0.33  0.57  1.06 1.54  2.67     3104     2460    1
#> d_ind            1.14 0.56  0.04  0.78  1.14 1.50  2.24     1464     2043    1
#> omega           -0.07 0.91 -1.83 -0.66 -0.08 0.50  1.85     2019     2461    1
#> tau              0.87 0.20  0.56  0.73  0.84 0.98  1.35     1197     1179    1
#> tau_consistency  0.84 0.18  0.55  0.71  0.82 0.95  1.25     1335     2260    1
#> 
#> Residual deviance: 53.8 (on 50 data points)
#>                pD: 43.9
#>               DIC: 97.7
#> 
#> Bayesian p-value: 0.93
#> 
#> ------------------------- Node-split Individual counselling vs. No intervention ---- 
#> 
#>                 mean   sd  2.5%   25%  50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net           0.85 0.24  0.41  0.69 0.84 0.99  1.34     1199     1957    1
#> d_dir           0.88 0.26  0.40  0.70 0.87 1.04  1.42     1983     2322    1
#> d_ind           0.59 0.69 -0.74  0.13 0.57 1.02  1.96     1574     1640    1
#> omega           0.29 0.72 -1.12 -0.17 0.28 0.77  1.70     1598     2070    1
#> tau             0.86 0.20  0.56  0.72 0.84 0.97  1.34     1285     1977    1
#> tau_consistency 0.84 0.18  0.55  0.71 0.82 0.95  1.25     1335     2260    1
#> 
#> Residual deviance: 54 (on 50 data points)
#>                pD: 44
#>               DIC: 98
#> 
#> Bayesian p-value: 0.65
#> 
#> -------------------------------------- Node-split Self-help vs. No intervention ---- 
#> 
#>                  mean   sd  2.5%   25%   50%  75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net            0.49 0.40 -0.27  0.22  0.48 0.76  1.34     1987     2363 1.00
#> d_dir            0.32 0.55 -0.81 -0.03  0.33 0.69  1.40     3691     2732 1.00
#> d_ind            0.70 0.62 -0.54  0.28  0.70 1.10  1.92     2162     2598 1.00
#> omega           -0.37 0.83 -2.05 -0.91 -0.35 0.18  1.24     2391     2271 1.00
#> tau              0.86 0.20  0.56  0.72  0.83 0.97  1.34     1233     2120 1.01
#> tau_consistency  0.84 0.18  0.55  0.71  0.82 0.95  1.25     1335     2260 1.00
#> 
#> Residual deviance: 54.2 (on 50 data points)
#>                pD: 44.5
#>               DIC: 98.7
#> 
#> Bayesian p-value: 0.65
#> 
#> ----------------------- Node-split Individual counselling vs. Group counselling ---- 
#> 
#>                  mean   sd  2.5%   25%   50%   75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net           -0.27 0.41 -1.05 -0.55 -0.27  0.00  0.53     3123     2385    1
#> d_dir           -0.12 0.49 -1.10 -0.42 -0.12  0.20  0.89     3380     3073    1
#> d_ind           -0.53 0.64 -1.81 -0.94 -0.52 -0.12  0.74     1113     1443    1
#> omega            0.42 0.70 -0.95 -0.03  0.41  0.87  1.81     1246     1725    1
#> tau              0.87 0.20  0.55  0.73  0.85  0.99  1.33     1023     1610    1
#> tau_consistency  0.84 0.18  0.55  0.71  0.82  0.95  1.25     1335     2260    1
#> 
#> Residual deviance: 53.8 (on 50 data points)
#>                pD: 44.4
#>               DIC: 98.2
#> 
#> Bayesian p-value: 0.53
#> 
#> ------------------------------------ Node-split Self-help vs. Group counselling ---- 
#> 
#>                  mean   sd  2.5%   25%   50%   75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net           -0.63 0.50 -1.62 -0.95 -0.62 -0.30  0.34     3109     2888    1
#> d_dir           -0.62 0.68 -1.97 -1.05 -0.60 -0.18  0.72     3889     2936    1
#> d_ind           -0.63 0.68 -1.98 -1.08 -0.63 -0.18  0.71     2097     2354    1
#> omega            0.01 0.90 -1.80 -0.56  0.01  0.59  1.82     2212     2436    1
#> tau              0.88 0.19  0.57  0.74  0.86  0.99  1.34     1446     2530    1
#> tau_consistency  0.84 0.18  0.55  0.71  0.82  0.95  1.25     1335     2260    1
#> 
#> Residual deviance: 53.8 (on 50 data points)
#>                pD: 44.3
#>               DIC: 98
#> 
#> Bayesian p-value: 0.99
#> 
#> ------------------------------- Node-split Self-help vs. Individual counselling ---- 
#> 
#>                  mean   sd  2.5%   25%   50%   75% 97.5% Bulk_ESS Tail_ESS Rhat
#> d_net           -0.36 0.41 -1.15 -0.62 -0.35 -0.09  0.46     2369     2496 1.00
#> d_dir            0.08 0.65 -1.21 -0.34  0.09  0.49  1.33     3636     2885 1.00
#> d_ind           -0.61 0.53 -1.69 -0.94 -0.60 -0.26  0.38     1826     2268 1.00
#> omega            0.69 0.82 -0.88  0.16  0.67  1.21  2.31     2230     2456 1.00
#> tau              0.86 0.19  0.56  0.72  0.83  0.97  1.28     1190     1885 1.01
#> tau_consistency  0.84 0.18  0.55  0.71  0.82  0.95  1.25     1335     2260 1.00
#> 
#> Residual deviance: 53.7 (on 50 data points)
#>                pD: 44.1
#>               DIC: 97.8
#> 
#> Bayesian p-value: 0.39
# }

## Plaque psoriasis ML-NMR
# Set up plaque psoriasis network combining IPD and AgD
library(dplyr)
pso_ipd <- filter(plaque_psoriasis_ipd,
                  studyc %in% c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"))

pso_agd <- filter(plaque_psoriasis_agd,
                  studyc == "FIXTURE")

head(pso_ipd)
#>      studyc      trtc_long    trtc trtn pasi75 pasi90 pasi100 age  bmi pasi_w0
#> 1 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      0      0       0  34 32.2    18.2
#> 2 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      1      0       0  64 41.9    23.4
#> 3 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      1      1       0  42 26.2    12.8
#> 4 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      0      0       0  45 52.9    36.0
#> 5 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      1      0       0  67 22.9    20.9
#> 6 UNCOVER-1 Ixekizumab Q2W IXE_Q2W    2      1      1       1  57 22.4    18.2
#>    male bsa weight durnpso prevsys   psa
#> 1  TRUE  18   98.1     6.7    TRUE  TRUE
#> 2  TRUE  33  129.6    14.5   FALSE  TRUE
#> 3  TRUE  33   78.0    26.5    TRUE FALSE
#> 4 FALSE  50  139.9    25.0    TRUE  TRUE
#> 5 FALSE  35   54.2    11.9    TRUE FALSE
#> 6  TRUE  29   67.5    15.2    TRUE FALSE
head(pso_agd)
#>    studyc          trtc_long    trtc trtn pasi75_r pasi75_n pasi90_r pasi90_n
#> 1 FIXTURE         Etanercept     ETN    4      142      323       67      323
#> 2 FIXTURE            Placebo     PBO    1       16      324        5      324
#> 3 FIXTURE Secukinumab 150 mg SEC_150    5      219      327      137      327
#> 4 FIXTURE Secukinumab 300 mg SEC_300    6      249      323      175      323
#>   pasi100_r pasi100_n sample_size_w0 age_mean age_sd bmi_mean bmi_sd
#> 1        14       323            326     43.8   13.0     28.7    5.9
#> 2         0       324            326     44.1   12.6     27.9    6.1
#> 3        47       327            327     45.4   12.9     28.4    5.9
#> 4        78       323            327     44.5   13.2     28.4    6.4
#>   pasi_w0_mean pasi_w0_sd male bsa_mean bsa_sd weight_mean weight_sd
#> 1         23.2        9.8 71.2     33.6   18.0        84.6      20.5
#> 2         24.1       10.5 72.7     35.2   19.1        82.0      20.4
#> 3         23.7       10.5 72.2     34.5   19.4        83.6      20.8
#> 4         23.9        9.9 68.5     34.3   19.2        83.0      21.6
#>   durnpso_mean durnpso_sd prevsys  psa
#> 1         16.4       12.0    65.6 13.5
#> 2         16.6       11.6    62.6 15.0
#> 3         17.3       12.2    64.8 15.0
#> 4         15.8       12.3    63.0 15.3

pso_ipd <- pso_ipd %>%
  mutate(# Variable transformations
    bsa = bsa / 100,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    weight = weight / 10,
    durnpso = durnpso / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker"),
    # Check complete cases for covariates of interest
    complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  )

pso_agd <- pso_agd %>%
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    prevsys = prevsys / 100,
    psa = psa / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker")
  )

# Exclude small number of individuals with missing covariates
pso_ipd <- filter(pso_ipd, complete)

pso_net <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = pasi75,
          trt_class = trtclass),
  set_agd_arm(pso_agd,
              study = studyc,
              trt = trtc,
              r = pasi75_r,
              n = pasi75_n,
              trt_class = trtclass)
)

# Print network details
pso_net
#> A network with 3 IPD studies, and 1 AgD study (arm-based).
#> 
#> ------------------------------------------------------------------- IPD studies ---- 
#>  Study     Treatment arms                  
#>  UNCOVER-1 3: IXE_Q2W | IXE_Q4W | PBO      
#>  UNCOVER-2 4: ETN | IXE_Q2W | IXE_Q4W | PBO
#>  UNCOVER-3 4: ETN | IXE_Q2W | IXE_Q4W | PBO
#> 
#>  Outcome type: binary
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study   Treatment arms                  
#>  FIXTURE 4: PBO | ETN | SEC_150 | SEC_300
#> 
#>  Outcome type: count
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 6, in 3 classes
#> Total number of studies: 4
#> Reference treatment is: PBO
#> Network is connected

# Add integration points to the network
pso_net <- add_integration(pso_net,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
  weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
  psa = distr(qbern, prob = psa),
  n_int = 64)
#> Using weighted average correlation matrix computed from IPD studies.

# \donttest{
# Fitting a ML-NMR model.
# Specify a regression model to include effect modifier interactions for five
# covariates, along with main (prognostic) effects. We use a probit link and
# specify that the two-parameter Binomial approximation for the aggregate-level
# likelihood should be used. We set treatment-covariate interactions to be equal
# within each class. We narrow the possible range for random initial values with
# init_r = 0.1, since probit models in particular are often hard to initialise.
# Using the QR decomposition greatly improves sampling efficiency here, as is
# often the case for regression models.
pso_fit <- nma(pso_net, 
               trt_effects = "fixed",
               link = "probit",
               likelihood = "bernoulli2",
               regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
               class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10),
               init_r = 0.1,
               QR = TRUE)
#> Note: Setting "PBO" as the network reference treatment.
pso_fit
#> A fixed effects ML-NMR with a bernoulli2 likelihood (probit link).
#> Regression model: ~(durnpso + prevsys + bsa + weight + psa) * .trt.
#> Centred covariates at the following overall mean values:
#>   durnpso   prevsys       bsa    weight       psa 
#> 1.8134535 0.6450416 0.2909089 8.9369318 0.2147914 
#> Inference for Stan model: binomial_2par.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>                                         mean se_mean   sd     2.5%      25%
#> beta[durnpso]                           0.04    0.00 0.06    -0.08     0.00
#> beta[prevsys]                          -0.14    0.00 0.16    -0.45    -0.25
#> beta[bsa]                              -0.07    0.01 0.45    -1.00    -0.36
#> beta[weight]                            0.04    0.00 0.03    -0.02     0.02
#> beta[psa]                              -0.08    0.00 0.17    -0.42    -0.19
#> beta[durnpso:.trtclassTNFa blocker]    -0.03    0.00 0.07    -0.17    -0.08
#> beta[durnpso:.trtclassIL blocker]      -0.01    0.00 0.07    -0.14    -0.06
#> beta[prevsys:.trtclassTNFa blocker]     0.19    0.00 0.19    -0.17     0.06
#> beta[prevsys:.trtclassIL blocker]       0.06    0.00 0.17    -0.27    -0.06
#> beta[bsa:.trtclassTNFa blocker]         0.06    0.01 0.52    -0.93    -0.29
#> beta[bsa:.trtclassIL blocker]           0.29    0.01 0.50    -0.67    -0.05
#> beta[weight:.trtclassTNFa blocker]     -0.17    0.00 0.04    -0.24    -0.19
#> beta[weight:.trtclassIL blocker]       -0.10    0.00 0.03    -0.16    -0.12
#> beta[psa:.trtclassTNFa blocker]        -0.05    0.00 0.21    -0.45    -0.19
#> beta[psa:.trtclassIL blocker]           0.01    0.00 0.18    -0.35    -0.11
#> d[ETN]                                  1.55    0.00 0.08     1.40     1.50
#> d[IXE_Q2W]                              2.95    0.00 0.09     2.78     2.89
#> d[IXE_Q4W]                              2.54    0.00 0.08     2.38     2.48
#> d[SEC_150]                              2.14    0.00 0.11     1.92     2.07
#> d[SEC_300]                              2.45    0.00 0.12     2.22     2.37
#> lp__                                -1576.37    0.09 3.46 -1584.16 -1578.48
#>                                          50%      75%    97.5% n_eff Rhat
#> beta[durnpso]                           0.04     0.09     0.16  5196    1
#> beta[prevsys]                          -0.14    -0.03     0.16  4999    1
#> beta[bsa]                              -0.05     0.24     0.79  4920    1
#> beta[weight]                            0.04     0.06     0.10  5407    1
#> beta[psa]                              -0.08     0.03     0.24  5913    1
#> beta[durnpso:.trtclassTNFa blocker]    -0.03     0.02     0.12  5812    1
#> beta[durnpso:.trtclassIL blocker]      -0.01     0.03     0.12  6260    1
#> beta[prevsys:.trtclassTNFa blocker]     0.19     0.32     0.55  5327    1
#> beta[prevsys:.trtclassIL blocker]       0.06     0.18     0.40  5968    1
#> beta[bsa:.trtclassTNFa blocker]         0.05     0.39     1.13  5376    1
#> beta[bsa:.trtclassIL blocker]           0.27     0.62     1.33  5975    1
#> beta[weight:.trtclassTNFa blocker]     -0.17    -0.14    -0.10  5261    1
#> beta[weight:.trtclassIL blocker]       -0.10    -0.08    -0.04  6266    1
#> beta[psa:.trtclassTNFa blocker]        -0.05     0.09     0.36  5963    1
#> beta[psa:.trtclassIL blocker]           0.01     0.13     0.38  6671    1
#> d[ETN]                                  1.55     1.60     1.71  4573    1
#> d[IXE_Q2W]                              2.95     3.01     3.13  4767    1
#> d[IXE_Q4W]                              2.54     2.60     2.70  5233    1
#> d[SEC_150]                              2.14     2.22     2.37  5488    1
#> d[SEC_300]                              2.45     2.53     2.68  5719    1
#> lp__                                -1576.05 -1573.88 -1570.57  1619    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 17 08:17:44 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }

## Newly-diagnosed multiple myeloma NMA
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
#> d[Len]                 -0.54    0.00 0.05    -0.62    -0.57    -0.54    -0.50
#> d[Thal]                -0.11    0.00 0.09    -0.29    -0.17    -0.11    -0.05
#> lp__                -6230.00    0.06 2.46 -6235.50 -6231.49 -6229.69 -6228.20
#> shape[Attal2012]        1.30    0.00 0.06     1.18     1.26     1.30     1.34
#> shape[Jackson2019]      0.93    0.00 0.02     0.89     0.92     0.93     0.95
#> shape[McCarthy2012]     1.29    0.00 0.07     1.16     1.25     1.29     1.34
#> shape[Morgan2012]       0.88    0.00 0.03     0.82     0.86     0.88     0.90
#> shape[Palumbo2014]      1.01    0.00 0.07     0.88     0.97     1.01     1.06
#>                        97.5% n_eff Rhat
#> d[Len]                 -0.45  5785    1
#> d[Thal]                 0.05  5270    1
#> lp__                -6226.20  1928    1
#> shape[Attal2012]        1.42  4507    1
#> shape[Jackson2019]      0.98  4813    1
#> shape[McCarthy2012]     1.43  4602    1
#> shape[Morgan2012]       0.94  5126    1
#> shape[Palumbo2014]      1.16  4504    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Apr 17 08:21:10 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
