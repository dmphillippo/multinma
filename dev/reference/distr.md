# Specify a general marginal distribution

`distr()` is used within the function
[`add_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md)
to specify marginal distributions for the covariates, via a
corresponding inverse CDF. It is also used in
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md)
to specify a distribution for the baseline response (intercept) when
predicting absolute outcomes.

## Usage

``` r
distr(qfun, ...)
```

## Arguments

- qfun:

  an inverse CDF, either as a function name or a string

- ...:

  parameters of the distribution as arguments to `qfun`, these will be
  quoted and evaluated later in the context of the aggregate data
  sources

## Value

An object of class distr.

## Details

The function `qfun` should have a formal argument called `p`. This
restriction serves as a crude check for inverse CDFs (e.g. an error will
be given if `dnorm` is used instead of `qnorm`). If a user-written CDF
is supplied, it must have an argument `p` which takes a vector of
probabilities.

## See also

[`add_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md)
where `distr()` is used to specify marginal distributions for covariates
to integrate over, and
[`predict.stan_nma()`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md)
where `distr()` is used to specify a distribution on the baseline
response.

## Examples

``` r
## Specifying marginal distributions for integration

df <- data.frame(x1_mean = 2, x1_sd = 0.5, x2 = 0.8)

# Distribution parameters are evaluated in the context of the data frame
add_integration(df,
                x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                x2 = distr(qbern, prob = x2),
                cor = diag(2))
#> # A tibble: 1 × 5
#>   x1_mean x1_sd    x2 .int_x1    .int_x2   
#>     <dbl> <dbl> <dbl> <list>     <list>    
#> 1       2   0.5   0.8 <dbl [64]> <dbl [64]>
```
