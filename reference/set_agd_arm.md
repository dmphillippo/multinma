# Set up arm-based aggregate data

Set up a network containing arm-based aggregate data (AgD), such as
event counts or mean outcomes on each arm. Multiple data sources may be
combined once created using
[`combine_network()`](https://dmphillippo.github.io/multinma/reference/combine_network.md).

## Usage

``` r
set_agd_arm(
  data,
  study,
  trt,
  y = NULL,
  se = NULL,
  r = NULL,
  n = NULL,
  E = NULL,
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

- r:

  column of `data` specifying a binary or Binomial outcome count

- n:

  column of `data` specifying Binomial outcome numerator

- E:

  column of `data` specifying the total time at risk for Poisson
  outcomes

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

The `sample_size` argument is optional, but when specified:

- Enables automatic centering of predictors (`center = TRUE`) in
  [`nma()`](https://dmphillippo.github.io/multinma/reference/nma.md)
  when a regression model is given for a network combining IPD and AgD

- Enables production of study-specific relative effects, rank
  probabilities, etc. for studies in the network when a regression model
  is given

- Nodes in
  [`plot.nma_data()`](https://dmphillippo.github.io/multinma/reference/plot.nma_data.md)
  may be weighted by sample size

If a Binomial outcome is specified and `sample_size` is omitted, `n`
will be used as the sample size by default. If a Multinomial outcome is
specified and `sample_size` is omitted, the sample size will be
determined automatically from the supplied counts by default.

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


# Plot network
plot(smk_net)
```
