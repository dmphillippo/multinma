# Set up contrast-based aggregate data

Set up a network containing contrast-based aggregate data (AgD), i.e.
summaries of relative effects between treatments such as log Odds
Ratios. Multiple data sources may be combined once created using
[`combine_network()`](https://dmphillippo.github.io/multinma/dev/reference/combine_network.md).

## Usage

``` r
set_agd_contrast(
  data,
  study,
  trt,
  y = NULL,
  se = NULL,
  sample_size = NULL,
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

- y:

  column of `data` specifying a continuous outcome

- se:

  column of `data` specifying the standard error for a continuous
  outcome

- sample_size:

  column of `data` giving the sample size in each arm. Optional, see
  details.

- trt_ref:

  reference treatment for the network, as a single integer, string, or
  factor. If not specified, a reasonable well-connected default will be
  chosen (see details).

- trt_class:

  column of `data` specifying treatment classes, coded using integers,
  strings, or factors. By default, no classes are specified.

## Value

An object of class
[nma_data](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)

## Details

Each study should have a single reference/baseline treatment, against
which relative effects in the other arm(s) are given. For the reference
arm, include a data row with continuous outcome `y` equal to `NA`. If a
study has three or more arms (so two or more relative effects), set the
standard error `se` for the reference arm data row equal to the standard
error of the mean outcome on the reference arm (this determines the
covariance of the relative effects, when expressed as differences in
mean outcomes between arms).

All arguments specifying columns of `data` accept the following:

- A column name as a character string, e.g. `study = "studyc"`

- A bare column name, e.g. `study = studyc`

- [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
  style semantics for inline variable transformations, e.g.
  `study = paste(author, year)`

By default, `trt_ref = NULL` and a network reference treatment will be
chosen that attempts to maximise computational efficiency and stability.
If an alternative reference treatment is chosen and the model runs
slowly or has low effective sample size (ESS) this may be the cause -
try letting the default reference treatment be used instead. Regardless
of which treatment is used as the network reference at the model fitting
stage, results can be transformed afterwards: see the `trt_ref` argument
of
[`relative_effects()`](https://dmphillippo.github.io/multinma/dev/reference/relative_effects.md)
and
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md).

The `sample_size` argument is optional, but when specified:

- Enables automatic centering of predictors (`center = TRUE`) in
  [`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
  when a regression model is given for a network combining IPD and AgD

- Enables production of study-specific relative effects, rank
  probabilities, etc. for studies in the network when a regression model
  is given

- Nodes in
  [`plot.nma_data()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_data.md)
  may be weighted by sample size

## See also

[`set_ipd()`](https://dmphillippo.github.io/multinma/dev/reference/set_ipd.md)
for individual patient data,
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md)
for arm-based aggregate data, and
[`combine_network()`](https://dmphillippo.github.io/multinma/dev/reference/combine_network.md)
for combining several data sources in one network.

[`print.nma_data()`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_data.md)
for the print method displaying details of the network, and
[`plot.nma_data()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_data.md)
for network plots.

## Examples

``` r
# Set up network of Parkinson's contrast data
head(parkinsons)
#>   studyn trtn     y    se   n  diff se_diff
#> 1      1    1 -1.22 0.504  54    NA   0.504
#> 2      1    3 -1.53 0.439  95 -0.31   0.668
#> 3      2    1 -0.70 0.282 172    NA   0.282
#> 4      2    2 -2.40 0.258 173 -1.70   0.382
#> 5      3    1 -0.30 0.505  76    NA   0.505
#> 6      3    2 -2.60 0.510  71 -2.30   0.718

park_net <- set_agd_contrast(parkinsons,
                             study = studyn,
                             trt = trtn,
                             y = diff,
                             se = se_diff,
                             sample_size = n)

# Print details
park_net
#> A network with 7 AgD studies (contrast-based).
#> 
#> -------------------------------------------------- AgD studies (contrast-based) ---- 
#>  Study Treatment arms
#>  1     2: 1 | 3      
#>  2     2: 1 | 2      
#>  3     3: 4 | 1 | 2  
#>  4     2: 4 | 3      
#>  5     2: 4 | 3      
#>  6     2: 4 | 5      
#>  7     2: 4 | 5      
#> 
#>  Outcome type: continuous
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 5
#> Total number of studies: 7
#> Reference treatment is: 4
#> Network is connected

# Plot network
plot(park_net)
```
