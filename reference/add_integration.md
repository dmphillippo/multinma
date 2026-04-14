# Add numerical integration points to aggregate data

The `add_integration()` generic creates Quasi-Monte Carlo numerical
integration points using a Gaussian copula and Sobol' sequences, as
described in Phillippo et al. (2020) . Methods are available for
networks stored in `nma_data` objects, and for data frames. The function
`unnest_integration()` unnests integration points stored in a data
frame, to aid plotting or other exploration.

## Usage

``` r
add_integration(x, ...)

# Default S3 method
add_integration(x, ...)

# S3 method for class 'data.frame'
add_integration(
  x,
  ...,
  cor = NULL,
  cor_adjust = NULL,
  n_int = 64L,
  int_args = list()
)

# S3 method for class 'nma_data'
add_integration(
  x,
  ...,
  cor = NULL,
  cor_adjust = NULL,
  n_int = 64L,
  int_args = list()
)

unnest_integration(data)
```

## Arguments

- x:

  An `nma_data` object, as created by the `set_*()` functions or
  [`combine_network()`](https://dmphillippo.github.io/multinma/reference/combine_network.md),
  or data frame

- ...:

  Distributions for covariates, see "Details"

- cor:

  Correlation matrix to use for generating the integration points. By
  default, this takes a weighted correlation matrix from all IPD
  studies. Rows and columns should match the order of covariates
  specified in `...`.

- cor_adjust:

  Adjustment to apply to the correlation matrix given by `cor` (or
  computed from the IPD if `cor = NULL`) to obtain the Gaussian copula
  correlations, either `"spearman"`, `"pearson"`, or `"none"`, see
  "Details". The default when `cor = NULL` is `"spearman"`, otherwise
  the default is `"pearson"`.

- n_int:

  Number of integration points to generate, default 64. Powers of 2 are
  recommended, as these are expected to be particularly efficient for
  QMC integration.

- int_args:

  A named list of arguments to pass to
  [`sobol()`](https://rdrr.io/pkg/randtoolbox/man/quasiRNG.html)

- data:

  Data frame with nested integration points, stored in list columns as
  `.int_<variable name>`

## Value

For the `nma_data` method, an object of class
[nma_data](https://dmphillippo.github.io/multinma/reference/nma_data-class.md).
For the `data.frame` method, the input data frame is returned (as a
`tibble`) with an added column for each covariate (prefixed with
".int\_"), containing the numerical integration points nested as
length-`n_int` vectors within each row. For `unnest_integration()`, a
data frame with integration points unnested.

## Details

The arguments passed to `...` specify distributions for the covariates.
Argument names specify the name of the covariate, which should match a
covariate name in the IPD (if IPD are present). The required marginal
distribution is then specified using the function
[`distr()`](https://dmphillippo.github.io/multinma/reference/distr.md).

The argument `cor_adjust` specifies how the correlation matrix given by
`cor` (or computed from the IPD if `cor = NULL`) is adjusted to obtain
the correlation matrix for the Gaussian copula, using the formulae in
Xiao and Zhou (2018) .

- `cor_adjust = "spearman"` should be used when the correlations `cor`
  have been computed using Spearman's rank correlation. Correlations
  between continuous covariates will be reproduced exactly by the
  integration points. Correlations between discrete covariates will be
  reproduced approximately. This is the default when `cor = NULL` and
  correlations are calculated from the IPD studies.

- `cor_adjust = "pearson"` should be used when the correlations `cor`
  have been computed using Pearson's product-moment correlation.
  Correlations between Normal covariates will be reproduced exactly by
  the integration points, all others will be reproduced approximately.
  Correlations between discrete covariates will be reproduced
  approximately (and identically to `cor_adjust = "spearman"`). This is
  the default when `cor` is provided by the user, since
  [`cor()`](https://rdrr.io/r/stats/cor.html) defaults to
  `method = "pearson"` and Pearson correlations are most likely reported
  in published data. However, we recommend providing Spearman
  correlations (e.g. from `cor(., method = "spearman")`) and using
  `cor_adjust = "spearman"` where possible.

- `cor_adjust = "none"` allows the user to specify the correlation
  matrix for the Gaussian copula directly; no adjustment is applied.

- `cor_adjust = "legacy"` is also available, which reproduces exactly
  the behaviour from version 0.3.0 and earlier. This is similar to
  `cor_adjust = "none"`, but unadjusted Spearman correlations are used
  if `cor = NULL`.

When adding integration points to a network object the correlation
matrix used is stored in `$int_cor`, and the copula correlation matrix
and adjustment used are stored as attributes of `$int_cor`. If this
correlation matrix is passed again to `add_integration()` (e.g. to reuse
the correlations for an external target population) this will be
detected, and the correct setting for `cor_adjust` will automatically be
applied.

## References

Phillippo DM, Dias S, Ades AE, Belger M, Brnabic A, Schacht A, Saure D,
Kadziola Z, Welton NJ (2020). “Multilevel Network Meta-Regression for
population-adjusted treatment comparisons.” *Journal of the Royal
Statistical Society: Series A (Statistics in Society)*, **183**(3),
1189–1210. [doi:10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579)
.  
  
Xiao Q, Zhou S (2018). “Matching a correlation coefficient by a Gaussian
copula.” *Communications in Statistics - Theory and Methods*, **48**(7),
1728–1747.
[doi:10.1080/03610926.2018.1439962](https://doi.org/10.1080/03610926.2018.1439962)
.

## Examples

``` r
## Plaque psoriasis ML-NMR - network setup and adding integration points
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


## Adding integration points to a data frame, e.g. for prediction
# Define a data frame of covariate summaries
new_agd_int <- data.frame(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1)

# Adding integration points, using the weighted average correlation matrix
# computed for the plaque psoriasis network
new_agd_int <- add_integration(new_agd_int,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
  weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
  psa = distr(qbern, prob = psa),
  cor = pso_net$int_cor,
  n_int = 64)

# Here, since we reused the correlation matrix pso_net$int_cor from the
# network, the correct setting of cor_adjust = "spearman" is automatically
# applied

new_agd_int
#> # A tibble: 1 × 13
#>   bsa_mean bsa_sd prevsys   psa weight_mean weight_sd durnpso_mean durnpso_sd
#>      <dbl>  <dbl>   <dbl> <dbl>       <dbl>     <dbl>        <dbl>      <dbl>
#> 1      0.6    0.3     0.1   0.2          10         1            3          1
#> # ℹ 5 more variables: .int_durnpso <list>, .int_prevsys <list>,
#> #   .int_bsa <list>, .int_weight <list>, .int_psa <list>
```
