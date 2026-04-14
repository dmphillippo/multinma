# The Gamma distribution

We provide convenient extensions of the `[dpq]gamma` functions, which
allow the distribution to be specified in terms of its mean and standard
deviation, instead of shape and rate/scale.

## Usage

``` r
qgamma(
  p,
  shape,
  rate = 1,
  scale = 1/rate,
  lower.tail = TRUE,
  log.p = FALSE,
  mean,
  sd
)

dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE, mean, sd)

pgamma(
  q,
  shape,
  rate = 1,
  scale = 1/rate,
  lower.tail = TRUE,
  log.p = FALSE,
  mean,
  sd
)
```

## Arguments

- p:

  vector of probabilities

- shape, rate, scale, log, lower.tail, log.p:

  see [stats::GammaDist](https://rdrr.io/r/stats/GammaDist.html)

- mean, sd:

  mean and standard deviation, overriding `shape` and `rate` or `scale`
  if specified

- x, q:

  vector of quantiles

## Value

Numeric vector of length equal to the maximum of the lengths of the
input arguments.
