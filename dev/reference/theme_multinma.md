# Plot theme for multinma plots

A simple `ggplot2` theme for plots in the `multinma` package.

## Usage

``` r
theme_multinma(...)
```

## Arguments

- ...:

  Arguments passed to
  [`ggplot2::theme_light()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)

## Value

A `ggplot2` theme

## See also

[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html),
[`ggplot2::theme_set()`](https://ggplot2.tidyverse.org/reference/get_theme.html)

## Examples

``` r
library(ggplot2)
theme_set(theme_multinma())
```
