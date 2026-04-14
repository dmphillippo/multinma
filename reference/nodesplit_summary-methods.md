# Methods for `nodesplit_summary` objects

The [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
and `as.tibble()` methods return the node-splitting summaries in a data
frame or tibble.

## Usage

``` r
# S3 method for class 'nodesplit_summary'
print(x, ..., digits = 2)

# S3 method for class 'nodesplit_summary'
as_tibble(x, ..., nest = FALSE)

as.tibble.nodesplit_summary(x, ..., nest = FALSE)

# S3 method for class 'nodesplit_summary'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `nodesplit_summary` object

- ...:

  Additional arguments passed on to other methods

- digits:

  Integer number of digits to display

- nest:

  Whether to return a nested tibble, with the full
  [nma_summary](https://dmphillippo.github.io/multinma/reference/nma_summary-class.md)
  and
  [nma_dic](https://dmphillippo.github.io/multinma/reference/nma_dic-class.md)
  objects, or to unnest their summaries, default `FALSE`

## Value

A `data.frame` for
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), a
`tbl_df` for `as.tibble()` and
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html).

The [`print()`](https://rdrr.io/r/base/print.html) method returns `x`
invisibly.

## See also

[`plot.nodesplit_summary()`](https://dmphillippo.github.io/multinma/reference/plot.nodesplit_summary.md)
