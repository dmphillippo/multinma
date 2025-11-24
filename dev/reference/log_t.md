# Log Student's t distribution

Density, distribution, and quantile function for the log t distribution,
whose logarithm has degrees of freedom `df`, mean `location`, and
standard deviation `scale`.

## Usage

``` r
dlogt(x, df, location = 0, scale = 1)

plogt(q, df, location = 0, scale = 1)

qlogt(p, df, location = 0, scale = 1)
```

## Arguments

- x, q:

  Vector of quantiles

- df:

  Degrees of freedom, greater than zero

- location:

  Location parameter

- scale:

  Scale parameter, greater than zero

- p:

  Vector of probabilities

## Value

`dlogt()` gives the density, `plogt()` gives the distribution function,
`qlogt()` gives the quantile function.

## Details

If \\\log(Y) \sim t\_\nu(\mu, \sigma^2)\\, then \\Y\\ has a log t
distribution with `location` \\\mu\\, `scale` \\\sigma\\, and `df`
\\\nu\\.

The mean and all higher moments of the log t distribution are undefined
or infinite.

If `df = 1` then the distribution is a log Cauchy distribution. As `df`
tends to infinity, this approaches a log Normal distribution.
