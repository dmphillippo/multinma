# Mean off-time reduction in Parkison's disease

Data frame containing the mean off-time reduction in patients given
dopamine agonists as adjunct therapy in Parkinson's disease, from 7
trials comparing four active drugs and placebo (Dias et al. 2011) .

## Usage

``` r
parkinsons
```

## Format

A data frame with 15 rows and 7 variables:

- studyn:

  numeric study ID

- trtn:

  numeric treatment code (placebo = 1)

- y:

  mean off-time reduction

- se:

  standard error

- n:

  sample size

- diff:

  mean difference vs. treatment in reference arm

- se_diff:

  standard error of mean difference, see details

## Details

This dataset may be analysed using either an arm-based likelihood using
`y` and `se`, or a contrast-based likelihood using `diff` and `se_diff`
(or a combination of the two across different studies).

The contrast-based data is formatted as described in
[`set_agd_contrast()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_contrast.md).
That is, for the chosen reference arm in each study, the mean difference
`diff` is set to `NA`, and `se_diff` is set to the standard error `se`
of the outcome on the reference arm.

## References

Dias S, Welton NJ, Sutton AJ, Ades AE (2011). “NICE DSU Technical
Support Document 2: A generalised linear modelling framework for
pair-wise and network meta-analysis of randomised controlled trials.”
National Institute for Health and Care Excellence.
<https://sheffield.ac.uk/nice-dsu>.
