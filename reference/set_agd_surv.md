# Set up aggregate survival data

Set up a network containing aggregate survival data (AgD) in the form of
event/censoring times (e.g. reconstructed from digitized Kaplan-Meier
curves) and covariate summary statistics from each study. Multiple data
sources may be combined once created using
[`combine_network()`](https://dmphillippo.github.io/multinma/reference/combine_network.md).

## Usage

``` r
set_agd_surv(
  data,
  study,
  trt,
  Surv,
  covariates = NULL,
  trt_ref = NULL,
  trt_class = NULL
)
```

## Arguments

- data:

  a data frame

- study:

  column of `data` specifying the studies, coded using integers,
  strings, or factors

- trt:

  column of `data` specifying treatments, coded using integers, strings,
  or factors

- Surv:

  column of `data` specifying a survival or time-to-event outcome, using
  the [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) function.
  Right/left/interval censoring and left truncation (delayed entry) are
  supported.

- covariates:

  data frame of covariate summary statistics for each study or study
  arm, with corresponding `study` and `trt` columns to match to those in
  `data`

- trt_ref:

  reference treatment for the network, as a single integer, string, or
  factor. If not specified, a reasonable well-connected default will be
  chosen (see details).

- trt_class:

  column of `data` specifying treatment classes, coded using integers,
  strings, or factors. By default, no classes are specified.

## Value

An object of class
[nma_data](https://dmphillippo.github.io/multinma/reference/nma_data-class.md)

## Details

By default, `trt_ref = NULL` and a network reference treatment will be
chosen that attempts to maximise computational efficiency and stability.
If an alternative reference treatment is chosen and the model runs
slowly or has low effective sample size (ESS) this may be the cause -
try letting the default reference treatment be used instead. Regardless
of which treatment is used as the network reference at the model fitting
stage, results can be transformed afterwards: see the `trt_ref` argument
of
[`relative_effects()`](https://dmphillippo.github.io/multinma/reference/relative_effects.md)
and
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/reference/predict.stan_nma.md).

All arguments specifying columns of `data` accept the following:

- A column name as a character string, e.g. `study = "studyc"`

- A bare column name, e.g. `study = studyc`

- [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
  style semantics for inline variable transformations, e.g.
  `study = paste(author, year)`

## See also

[`set_ipd()`](https://dmphillippo.github.io/multinma/reference/set_ipd.md)
for individual patient data,
[`set_agd_contrast()`](https://dmphillippo.github.io/multinma/reference/set_agd_contrast.md)
for contrast-based aggregate data, and
[`combine_network()`](https://dmphillippo.github.io/multinma/reference/combine_network.md)
for combining several data sources in one network.

[`print.nma_data()`](https://dmphillippo.github.io/multinma/reference/print.nma_data.md)
for the print method displaying details of the network, and
[`plot.nma_data()`](https://dmphillippo.github.io/multinma/reference/plot.nma_data.md)
for network plots.

## Examples

``` r
## Newly diagnosed multiple myeloma

head(ndmm_agd)  # Reconstructed Kaplan-Meier data
#>        study     studyf trt trtf eventtime status
#> 1 Morgan2012 Morgan2012 Pbo  Pbo  18.72575      1
#> 2 Morgan2012 Morgan2012 Pbo  Pbo  63.36000      0
#> 3 Morgan2012 Morgan2012 Pbo  Pbo  34.35726      1
#> 4 Morgan2012 Morgan2012 Pbo  Pbo  10.77826      1
#> 5 Morgan2012 Morgan2012 Pbo  Pbo  63.36000      0
#> 6 Morgan2012 Morgan2012 Pbo  Pbo  14.52966      1
ndmm_agd_covs   # Summary covariate information on each arm
#>         study      studyf  trt trtf sample_size  age_min age_iqr_l age_median
#> 1 Jackson2019 Jackson2019  Len  Len        1137 17.28246  59.13164   65.76766
#> 2 Jackson2019 Jackson2019  Pbo  Pbo         864 21.18572  58.30991   65.47402
#> 3  Morgan2012  Morgan2012  Pbo  Pbo         410 33.88979  58.05696   64.15999
#> 4  Morgan2012  Morgan2012 Thal Thal         408 38.45127  59.30022   65.48736
#>   age_iqr_h  age_max age_mean   age_sd iss_stage3 response_cr_vgpr      male
#> 1  72.00756 85.76095 65.16867 8.936962  0.2480211        0.8258575 0.6165347
#> 2  71.80261 86.23080 64.62894 9.399272  0.1921296        0.8310185 0.6215278
#> 3  70.44791 84.79372 63.92360 9.006311  0.3634146        0.7170732 0.6195122
#> 4  71.73597 84.69365 65.59387 8.384686  0.3186275        0.7450980 0.6151961

set_agd_surv(ndmm_agd,
             study = studyf,
             trt = trtf,
             Surv = Surv(eventtime, status),
             covariates = ndmm_agd_covs)
#> A network with 2 AgD studies (arm-based).
#> 
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study       Treatment arms
#>  Jackson2019 2: Pbo | Len  
#>  Morgan2012  2: Pbo | Thal 
#> 
#>  Outcome type: survival
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 3
#> Total number of studies: 2
#> Reference treatment is: Pbo
#> Network is connected
```
