# The Bernoulli Distribution

The density function `dbern()`, distribution function `pbern()`, and
quantile function `qbern()` for a Bernoulli distribution, with success
probability `prob`. These are equivalent to `dbinom(p, 1, prob)`,
`pbinom(p, 1, prob)` and `qbinom(p, 1, prob)`.

## Usage

``` r
qbern(p, prob, lower.tail = TRUE, log.p = FALSE)

pbern(q, prob, lower.tail = TRUE, log.p = FALSE)

dbern(x, prob, log = FALSE)
```

## Arguments

- p:

  vector of probabilities

- prob:

  probability of success

- lower.tail, log.p, log:

  see [stats::Binomial](https://rdrr.io/r/stats/Binomial.html)

- x, q:

  vector of quantiles

## Value

Numeric vector of length equal to the maximum of the lengths of the
input arguments.
