# HTA Plaque Psoriasis

Data frame containing the results of 16 trials comparing 8 treatments
for moderate-to-severe plaque psoriasis from an HTA report (Woolacott et
al. 2006) , analysed in TSD2 (Dias et al. 2011) . Outcomes are
success/failure to achieve 50%, 75%, or 90% reduction in symptoms on the
Psoriasis Area and Severity Index (PASI) scale. Some studies report all
three ordered outcomes, others only one or two. The latter are coded as
missing values (see details).

## Usage

``` r
hta_psoriasis
```

## Format

A data frame with 36 rows and 9 variables:

- studyn:

  numeric study ID

- studyc:

  study name

- year:

  year of publication

- trtn:

  numeric treatment code

- trtc:

  treatment name

- sample_size:

  sample size in each arm

- PASI50, PASI75, PASI90:

  ordered multinomial outcome counts (exclusive, see details)

## Details

Outcome counts are "exclusive"; that is, for a study reporting all
outcomes, the counts represent the categories 50 \< PASI \< 75, 75 \<
PASI \< 90, and 90 \< PASI \< 100, and are named by the lower end of the
interval. (As opposed to "inclusive" counts, which would represent the
overlapping categories PASI \> 50, PASI \> 70, and PASI \> 90.) The
count for the fourth category (the lowest), 0 \< PASI \< 50, is equal to
`sample_size - PASI50 - PASI75 - PASI90`.

Missing values are used where studies only report a subset of the
outcomes. For a study reporting only two outcomes, say 50 and 75, the
counts represent 50 \< PASI \< 75 and 75 \< PASI \< 100. For a study
reporting only one outcome, say PASI 75, the count represents 75 \< PASI
\< 100.

## References

Dias S, Welton NJ, Sutton AJ, Ades AE (2011). “NICE DSU Technical
Support Document 2: A generalised linear modelling framework for
pair-wise and network meta-analysis of randomised controlled trials.”
National Institute for Health and Care Excellence.
<https://www.sheffield.ac.uk/nice-dsu>.  
  
Woolacott N, Hawkins N, Mason A, Kainth A, Khadjesari Z, Bravo Vergel Y,
Misso K, Light K, Chalmers R, Sculpher M, Riemsma R (2006). “Etanercept
and efalizumab for the treatment of psoriasis: a systematic review.”
*Health Technology Assessment*, **10**(46).
[doi:10.3310/hta10460](https://doi.org/10.3310/hta10460) .
