# Generalised Student's t distribution (with location and scale)

Density, distribution, and quantile function for the generalised t
distribution with degrees of freedom `df`, shifted by `location` and
scaled by `scale`.

## Usage

``` r
dgent(x, df, location = 0, scale = 1)

pgent(q, df, location = 0, scale = 1)

qgent(p, df, location = 0, scale = 1)
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

`dgent()` gives the density, `pgent()` gives the distribution function,
`qgent()` gives the quantile function.
