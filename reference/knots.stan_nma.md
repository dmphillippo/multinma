# Knot locations for a fitted model

Obtain the knot locations from a fitted M-spline or piecewise
exponential model.

## Usage

``` r
# S3 method for class 'stan_nma'
knots(Fn, type = c("all", "internal", "boundary"), ...)
```

## Arguments

- Fn:

  A fitted
  [stan_nma](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
  object

- type:

  String, indicating whether to return all knots (`"all"`, the default),
  or only the internal knots (`"internal"`) or boundary knots
  (`"boundary"`).

- ...:

  Other arguments, passed on to
  [`splines2::knots()`](https://wwenjie.org/splines2/reference/knots.html)

## Value

A list of vectors of knot locations, for each study in the network.

## See also

[`make_knots()`](https://dmphillippo.github.io/multinma/reference/make_knots.md)
for constructing knots for a network object.
