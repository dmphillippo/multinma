# Methods for `nma_summary` objects

The [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
and `as.tibble()` methods return the posterior summary statistics in a
data frame or tibble. The
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) method returns a
matrix of posterior draws. The
[`as.array()`](https://rdrr.io/r/base/array.html) method returns a 3D
array \[Iteration, Chain, Parameter\] of posterior draws (as class
[mcmc_array](https://dmphillippo.github.io/multinma/reference/mcmc_array-class.md)).

## Usage

``` r
# S3 method for class 'nma_summary'
print(x, ..., digits = 2, pars, include = TRUE)

# S3 method for class 'nma_summary'
as.data.frame(x, ...)

# S3 method for class 'nma_summary'
as.tibble(x, ...)

# S3 method for class 'nma_summary'
as_tibble(x, ...)

# S3 method for class 'nma_summary'
as.array(x, ...)

# S3 method for class 'nma_summary'
as.matrix(x, ...)

# S3 method for class 'nma_rank_probs'
as.array(x, ...)

# S3 method for class 'nma_rank_probs'
as.matrix(x, ...)
```

## Arguments

- x:

  A `nma_summary` object

- ...:

  Additional arguments passed on to other methods

- digits:

  Integer number of digits to display

- pars:

  Character vector of parameters to display in the printed summary

- include:

  Logical, are parameters named in `pars` included (`TRUE`) or excluded
  (`FALSE`)

## Value

A `data.frame` for
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), a
`tbl_df` for `as.tibble()` and
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
a `matrix` for [`as.matrix()`](https://rdrr.io/r/base/matrix.html), and
an `mcmc_array` for [`as.array()`](https://rdrr.io/r/base/array.html).

The [`print()`](https://rdrr.io/r/base/print.html) method returns `x`
invisibly.

## See also

[`plot.nma_summary()`](https://dmphillippo.github.io/multinma/reference/plot.nma_summary.md)
