# Plaque psoriasis data

Two data frames, `plaque_psoriasis_ipd` and `plaque_psoriasis_agd`,
containing (simulated) individual patient data from four studies and
aggregate data from five studies (Phillippo 2019) . Outcomes are binary
success/failure to achieve 75%, 90%, or 100% reduction in symptoms on
the Psoriasis Area and Severity Index (PASI) scale.

## Usage

``` r
plaque_psoriasis_ipd

plaque_psoriasis_agd
```

## Format

The individual patient data are contained in a data frame
`plaque_psoriasis_ipd` with 4118 rows, one per individual, and 16
variables:

- studyc:

  study name

- trtc_long:

  treatment name (long format)

- trtc:

  treatment name

- trtn:

  numeric treatment code

- pasi75:

  binary PASI 75 outcome

- pasi90:

  binary PASI 90 outcome

- pasi100:

  binary PASI 100 outcome

- age:

  age (years)

- bmi:

  body mass index (BMI)

- pasi_w0:

  PASI score at week 0

- male:

  male sex (TRUE or FALSE)

- bsa:

  body surface area (percent)

- weight:

  weight (kilograms)

- durnpso:

  duration of psoriasis (years)

- prevsys:

  previous systemic treatment (TRUE or FALSE)

- psa:

  psoriatic arthritis (TRUE or FALSE)

The aggregate data are contained in a data frame `plaque_psoriasis_agd`
with 15 rows, one per study arm, and 26 variables:

- studyc:

  study name

- trtc_long:

  treatment name (long format)

- trtc:

  treatment name

- trtn:

  numeric treatment code

- pasi75_r, pasi75_n:

  PASI 75 outcome count and denominator

- pasi90_r, pasi90_n:

  PASI 75 outcome count and denominator

- pasi100_r, pasi100_n:

  PASI 75 outcome count and denominator

- sample_size_w0:

  sample size at week zero

- age_mean, age_sd:

  mean and standard deviation of age (years)

- bmi_mean, bmi_sd:

  mean and standard deviation of BMI

- pasi_w0_mean, pasi_w0_sd:

  mean and standard deviation of PASI score at week 0

- male:

  percentage of males

- bsa_mean, bsa_sd:

  mean and standard deviation of body surface area (percent)

- weight_mean, weight_sd:

  mean and standard deviation of weight (kilograms)

- durnpso_mean, durnpso_sd:

  mean and standard deviation of duration of psoriasis (years)

- prevsys:

  percentage of individuals with previous systemic treatment

- psa:

  percentage of individuals with psoriatic arthritis

An object of class `data.frame` with 15 rows and 26 columns.

## References

Phillippo DM (2019). *Calibration of Treatment Effects in Network
Meta-Analysis using Individual Patient Data*. Ph.D. thesis, University
of Bristol. Available from <https://research-information.bris.ac.uk/>.
