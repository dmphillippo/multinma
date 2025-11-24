# Summary of prior distributions

Print a summary of prior distribution details.

## Usage

``` r
# S3 method for class 'nma_prior'
summary(object, ..., probs = c(0.5, 0.95), digits = 2, trunc = NULL)
```

## Arguments

- object:

  Prior distribution as a `nma_prior` object

- ...:

  Additional arguments, not used

- probs:

  Numeric vector of probabilities to calculate prior intervals

- digits:

  Number of digits to display

- trunc:

  Optional numeric vector of length 2, giving the truncation limits of
  the prior distribution. Useful if a real-valued prior is assigned to a
  positive-valued parameter, then `trunc = c(0, Inf)` will give the
  correct prior intervals. By default, truncation is not used.

## Value

A data frame is returned invisibly, giving the prior intervals

## Examples

``` r
summary(normal(location = 0, scale = 1))
#> A Normal prior distribution: location = 0, scale = 1.
#> 50% of the prior density lies between -0.67 and 0.67.
#> 95% of the prior density lies between -1.96 and 1.96.
summary(half_normal(scale = 1))
#> A half-Normal prior distribution: scale = 1.
#> 50% of the prior density lies between 0 and 0.67.
#> 95% of the prior density lies between 0 and 1.96.
summary(log_normal(location = -3.93, scale = 1.51))
#> A log-Normal prior distribution: location = -3.93, scale = 1.51.
#> 50% of the prior density lies between 0.01 and 0.05.
#> 95% of the prior density lies between 0 and 0.38.

# Truncation limits may be set, for example to restrict a prior to positive values
summary(normal(location = 0.5, scale = 1), trunc = c(0, Inf))
#> A Normal prior distribution: location = 0.5, scale = 1.
#> 50% of the prior density lies between 0.45 and 1.44.
#> 95% of the prior density lies between 0.05 and 2.61.
```
