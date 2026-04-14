# Softmax transform

The softmax transform is a multivariate generalisation of the logit
transform. `softmax()` maps a vector of \\K-1\\ values on the real line
to a \\K\\-simplex (i.e. values between 0 and 1, that sum to 1).
`inv_softmax()` provides the inverse transform, mapping a \\K\\-simplex
vector to a vector of \\K-1\\ real values.

## Usage

``` r
softmax(x)

inv_softmax(p)
```

## Arguments

- x:

  \\K-1\\ vector of reals

- p:

  \\K\\ vector simplex

## Value

`softmax()` returns a vector of length \\K\\ that is a simplex.
`inv_softmax()` returns a vector of reals of length \\K-1\\.

## Examples

``` r
x <- c(-1, 3, -0.5, 2)
(p <- softmax(x))
#> [1] 0.03395701 0.01249208 0.68204471 0.02059597 0.25091023
sum(p)
#> [1] 1
inv_softmax(p)
#> [1] -1.0  3.0 -0.5  2.0
```
