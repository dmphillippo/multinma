# Newly diagnosed multiple myeloma

Three data frames, `ndmm_ipd`, `ndmm_agd`, and `ndmm_agd_covs`
containing (simulated) individual patient data (IPD) from three studies
and aggregate data (AgD) from two studies on newly diagnosed multiple
myeloma. The outcome of interest is progression-free survival after
autologous stem cell transplant. The IPD studies in `ndmm_ipd` provide
event/censoring times and covariate values for each individual. The AgD
studies provide reconstructed event/censoring times from digitized
Kaplan-Meier curves in `ndmm_agd` and covariate summaries in
`ndmm_agd_covs`, obtained from published trial reports. The data are
constructed to resemble those used by Leahy and Walsh (2019) .

## Usage

``` r
ndmm_ipd

ndmm_agd

ndmm_agd_covs
```

## Format

The individual patient data are contained in a data frame `ndmm_ipd`
with 1325 rows, one per individual, and 10 variables:

- study, studyf:

  study name

- trt, trtf:

  treatment name

- eventtime:

  event/censoring time

- status:

  censoring indicator (0 = censored, 1 = event)

- age:

  age (years)

- iss_stage3:

  ISS stage 3 (0 = no, 1 = yes)

- response_cr_vgpr:

  complete or very good partial response (0 = no, 1 = yes)

- male:

  male sex (0 = no, 1 = yes)

The reconstructed Kaplan-Meier data for the aggregate studies are
contained in a data frame `ndmm_agd` with 2819 rows and 6 variables:

- study, studyf:

  study name

- trt, trtf:

  treatment name

- eventtime:

  event/censoring time

- status:

  censoring indicator (0 = censored, 1 = event)

The covariate summaries extracted from published reportes for the
aggregate studies are contained in a data frame `ndmm_agd_covs` with 4
rows, one per study arm, and 15 columns:

- study, studyf:

  study name

- trt, trtf:

  treatment name

- sample_size:

  sample size in each arm

- age_min, age_iqr_l, age_median, age_iqr_h, age_max, age_mean, age_sd:

  summary statistics for age (years)

- iss_stage3:

  proportion of participants with ISS stage 3

- response_cr_vgpr:

  proportion of participants with complete or very good partial response

- male:

  proportion of male participants

## References

Leahy J, Walsh C (2019). “Assessing the impact of a matching-adjusted
indirect comparison in a Bayesian network meta-analysis.” *Research
Synthesis Methods*, **10**(4), 546–568.
[doi:10.1002/jrsm.1372](https://doi.org/10.1002/jrsm.1372) .
