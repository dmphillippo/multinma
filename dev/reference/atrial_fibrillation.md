# Stroke prevention in atrial fibrillation patients

Data frame containing the results of 26 trials comparing 17 treatments
in 4 classes for the prevention of stroke in patients with atrial
fibrillation (Cooper et al. 2009) . The data are the corrected versions
given by van Valkenhoef and Kuiper (2016) .

## Usage

``` r
atrial_fibrillation
```

## Format

A data frame with 63 rows and 11 variables:

- studyc:

  study name

- studyn:

  numeric study ID

- trtc:

  treatment name

- trtn:

  numeric treatment code

- trt_class:

  treatment class

- r:

  number of events

- n:

  sample size

- E:

  person-years at risk

- stroke:

  proportion of individuals with prior stroke

- year:

  year of study publication

- followup:

  mean length of follow-up (years)

## References

Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ (2009). “Addressing
between-study heterogeneity and inconsistency in mixed treatment
comparisons: Application to stroke prevention treatments in individuals
with non-rheumatic atrial fibrillation.” *Statistics in Medicine*,
**28**(14), 1861–1881.
[doi:10.1002/sim.3594](https://doi.org/10.1002/sim.3594) .  
  
van Valkenhoef G, Kuiper J (2016). *gemtc: Network Meta-Analysis Using
Bayesian Methods*. R package version 0.8-2,
<https://CRAN.R-project.org/package=gemtc>.
