# Multinomial outcome data

This function aids the specification of multinomial outcome data when
setting up a network with
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md)
or
[`set_ipd()`](https://dmphillippo.github.io/multinma/dev/reference/set_ipd.md).
It takes a set of columns (or, more generally, numeric vectors of the
same length) of outcome counts in each category, and binds these
together to produce a matrix.

## Usage

``` r
multi(..., inclusive = FALSE, type = c("ordered", "competing"))
```

## Arguments

- ...:

  Two or more numeric columns (or vectors) of category counts. Argument
  names (optional) will be used to label the categories.

- inclusive:

  Logical, are ordered category counts inclusive (`TRUE`) or exclusive
  (`FALSE`)? Default `FALSE`. Only used when `type = "ordered"`. See
  details.

- type:

  String, indicating whether categories are `"ordered"` or
  `"competing"`. Currently only ordered categorical outcomes are
  supported by the modelling functions in this package.

## Value

A matrix of (exclusive) category counts

## Details

When specifying ordered categorical counts, these can either be given as
*exclusive* counts (`inclusive = FALSE`, the default) where individuals
are only counted in the highest category they achieve, or *inclusive*
counts (`inclusive = TRUE`) where individuals are counted in every
category up to and including the highest category achieved. (Competing
outcomes, by nature, are always specified as exclusive counts.)

`NA` values can be used to indicate categories/cutpoints that were not
measured.

## Examples

``` r
# These two data sets specify the same ordered categorical data for outcomes
# r0 < r1 < r2, but the first uses the "inclusive" format and the second the
# "exclusive" format.
df_inclusive <- tibble::tribble(~r0, ~r1, ~r2,
                                1, 1, 1,
                                5, 4, 1,
                                5, 2, 2,
                                10, 5, 0,
                                5, 5, 0,
                                7, NA, 6,   # Achieved r2 or not (no r1)
                                10, 4, NA)  # Achieved r1 or not (no r2)

df_exclusive <- tibble::tribble(~r0, ~r1, ~r2,
                                0, 0, 1,
                                1, 3, 1,
                                3, 0, 2,
                                5, 5, 0,
                                0, 5, 0,
                                1, NA, 6,   # Achieved r2 or not (no r1)
                                6, 4, NA)   # Achieved r1 or not (no r2)

(r_inclusive <- with(df_inclusive, multi(r0, r1, r2, inclusive = TRUE)))
#>      r0 r1 r2
#> [1,]  0  0  1
#> [2,]  1  3  1
#> [3,]  3  0  2
#> [4,]  5  5  0
#> [5,]  0  5  0
#> [6,]  1 NA  6
#> [7,]  6  4 NA
#> attr(,"class")
#> [1] "multi_ordered" "matrix"        "array"        
(r_exclusive <- with(df_exclusive, multi(r0, r1, r2, inclusive = FALSE)))
#>      r0 r1 r2
#> [1,]  0  0  1
#> [2,]  1  3  1
#> [3,]  3  0  2
#> [4,]  5  5  0
#> [5,]  0  5  0
#> [6,]  1 NA  6
#> [7,]  6  4 NA
#> attr(,"class")
#> [1] "multi_ordered" "matrix"        "array"        

# Counts are always stored in exclusive format
stopifnot(isTRUE(all.equal(r_inclusive, r_exclusive)))


## HTA Plaque Psoriasis
library(dplyr)

# Ordered outcomes here are given as "exclusive" counts
head(hta_psoriasis)
#>   studyn   studyc year trtn             trtc sample_size PASI50 PASI75 PASI90
#> 1      1  Elewski 2004    1  Supportive care         193     12      5      1
#> 2      1  Elewski 2004    2 Etanercept 25 mg         196     59     46     21
#> 3      1  Elewski 2004    3 Etanercept 50 mg         194     54     56     40
#> 4      2 Gottlieb 2003    1  Supportive care          55      5      1      0
#> 5      2 Gottlieb 2003    2 Etanercept 25 mg          57     23     11      6
#> 6      3  Lebwohl 2003    1  Supportive care         122     13      5      1

# Calculate lowest category count (failure to achieve PASI 50)
pso_dat <- hta_psoriasis %>%
  mutate(`PASI<50` = sample_size - rowSums(cbind(PASI50, PASI75, PASI90), na.rm = TRUE))

# Set up network
pso_net <- set_agd_arm(pso_dat,
                       study = paste(studyc, year),
                       trt = trtc,
                       r = multi(`PASI<50`, PASI50, PASI75, PASI90,
                                 inclusive = FALSE,
                                 type = "ordered"))

pso_net
#> A network with 16 AgD studies (arm-based).
#> 
#> ------------------------------------------------------- AgD studies (arm-based) ---- 
#>  Study         Treatment arms                                          
#>  ACD2058g 2004 2: Supportive care | Efalizumab                         
#>  ACD2600g 2004 2: Supportive care | Efalizumab                         
#>  Altmeyer 1994 2: Supportive care | Fumaderm                           
#>  Chaudari 2001 2: Supportive care | Infliximab                         
#>  Elewski 2004  3: Supportive care | Etanercept 25 mg | Etanercept 50 mg
#>  Ellis 1991    3: Supportive care | Ciclosporin | Ciclosporin          
#>  Gordon 2003   2: Supportive care | Efalizumab                         
#>  Gottlieb 2003 2: Supportive care | Etanercept 25 mg                   
#>  Gottlieb 2004 3: Supportive care | Infliximab | Infliximab            
#>  Guenther 1991 2: Supportive care | Ciclosporin                        
#>  ... plus 6 more studies
#> 
#>  Outcome type: ordered (4 categories)
#> ------------------------------------------------------------------------------------
#> Total number of treatments: 8
#> Total number of studies: 16
#> Reference treatment is: Supportive care
#> Network is connected
```
