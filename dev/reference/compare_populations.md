# Compare population covariate distributions

Covariate distributions for studies in the network may be compared by
overlap effective sample size (via propensity scores) or by
(standardised) Euclidean distance between the covariate means. The
overlap effective sample size accounts for the full multivariate joint
distribution of the covariates, whereas the Euclidean distance between
the covariate means only accounts for the location and not the range of
variation or joint structure.

## Usage

``` r
compare_populations(
  network,
  covariates = NULL,
  method = c("propensity", "euclidean")
)

# S3 method for class 'pop_comp'
print(
  x,
  order = c("decreasing", "increasing"),
  simplify = TRUE,
  ...,
  digits = 2L,
  n = 10
)

# S3 method for class 'pop_comp'
plot(x, ...)
```

## Arguments

- network:

  An `nma_data` network object.

- covariates:

  Character vector of covariate names to compare on.

- method:

  Method to compare distributions, either `"propensity"` to calculate
  overlap effective sample size with propensity scores, or `"euclidean"`
  to calculate standardised Euclidean distance between means.

- x:

  A `pop_comp` object produced by `compare_populations()`

- order:

  String, the order in which to display the comparison summaries. Either
  `"decreasing"`, to list most similar studies first (the default), or
  `"increasing"` to list most dissimilar studies first.

- simplify:

  Logical, should the output be simplified to only show comparisons
  between subnetworks (`TRUE`, default), or between every study in the
  network (`FALSE`)? If a connected network was provided, then
  `simplify` is always `FALSE`.

- ...:

  Additional arguments (unused)

- digits:

  Number of digits to print in the summary, default 2

- n:

  Number of rows to show in the table of study summaries

## Value

A `pop_comp` object, containing a `summary` data frame and
`comparison_matrix` matrix of pairwise comparisons. When
`method = "propensity"`, a list `propensity_scores` of data frames of
fitted propensity scores will also be included.

## Aggregate data setup

When aggregate data are present in the network, some setup is necessary
to use this function. For `method = "propensity"`, the network must have
integration points present; specify these using
[`add_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md).
For `method = "euclidean"`, the network must either have integration
points present, or the mean and standard deviation of covariates can be
provided in the original input data when setting up the network. In the
latter case, the covariate means are assumed to be provided in the named
`<covariate>` column, and the corresponding standard deviation in a
column named `<covariate>_sd`. For example, if the covariate is `age`,
mean age should be provided in the `age` column, and the standard
devation of age in `age_sd` column.
