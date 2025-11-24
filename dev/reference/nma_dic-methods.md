# Methods for `nma_dic` objects

The [`print()`](https://rdrr.io/r/base/print.html) method prints details
of DIC model fit statistics, computed by the
[`dic()`](https://dmphillippo.github.io/multinma/dev/reference/dic.md)
function. The
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
and `as.tibble()` methods return the pointwise contributions to the DIC
and \$p_D\$ in a data frame or tibble. The
[`as.array()`](https://rdrr.io/r/base/array.html) and
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) methods returns a 3D
MCMC array (as class
[mcmc_array](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md))
or matrix of posterior draws of the residual deviances.

## Usage

``` r
# S3 method for class 'nma_dic'
print(x, digits = 1, ...)

# S3 method for class 'nma_dic'
as.data.frame(x, ...)

# S3 method for class 'nma_dic'
as.tibble(x, ...)

# S3 method for class 'nma_dic'
as_tibble(x, ...)

# S3 method for class 'nma_dic'
as.array(x, ...)

# S3 method for class 'nma_dic'
as.matrix(x, ...)
```

## Arguments

- x:

  An object of class
  [nma_dic](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-class.md)

- digits:

  Integer number of digits to display

- ...:

  Additional arguments passed on to other methods

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

[`dic()`](https://dmphillippo.github.io/multinma/dev/reference/dic.md),
[`plot.nma_dic()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_dic.md)
