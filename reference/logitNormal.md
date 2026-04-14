# The logit Normal distribution

Density, distribution function, and quantile function for the logit
Normal distribution. The location and scale parameters of the
distribution are `mu` and `sigma`, which are the mean and standard
deviation on the logit scale. For convenience, the distribution may also
be specified in terms of its mean and standard deviation, instead of its
logit-mean and logit-sd.

## Usage

``` r
dlogitnorm(x, mu = 0, sigma = 1, ..., mean, sd)

plogitnorm(q, mu = 0, sigma = 1, ..., mean, sd)

qlogitnorm(p, mu = 0, sigma = 1, ..., mean, sd)
```

## Arguments

- x, q:

  vector of quantiles

- mu, sigma:

  location and scale parameters, on the logit scale

- ...:

  additional arguments, passed to `[dpq]norm()`

- mean, sd:

  mean and standard deviation, overriding `mu` and `sigma` if specified

- p:

  vector of probabilities

## Value

Numeric vector of length equal to the maximum of the lengths of the
input arguments.
