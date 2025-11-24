# The nma_data class

The `nma_data` class contains the data for a NMA in a standard format,
created using the functions
[`set_ipd()`](https://dmphillippo.github.io/multinma/dev/reference/set_ipd.md),
[`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md),
[`set_agd_contrast()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_contrast.md),
[`set_agd_surv()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_surv.md),
or
[`combine_network()`](https://dmphillippo.github.io/multinma/dev/reference/combine_network.md).
The sub-class `mlnmr_data` is created by the function
[`add_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md),
and further contains numerical integration points for the aggregate
data.

## Details

Objects of class `nma_data` have the following components:

- `agd_arm`:

  data from studies with aggregate data (arm format)

- `agd_contrast`:

  data from studies with aggregate data (contrast format)

- `ipd`:

  data from studies with individual patient data

- `treatments`:

  treatment coding factor for entire network

- `classes`:

  treatment class coding factor (same length as `treatments` for entire
  network)

- `studies`:

  study coding factor for entire network

- `outcome`:

  outcome type for each data source, named list

The `agd_arm`, `agd_contrast`, and `ipd` components are tibbles with the
following columns:

- `.study`:

  study (as factor)

- `.trt`:

  treatment (as factor)

- `.trtclass`:

  treatment class (as factor), if specified

- `.y`:

  continuous outcome

- `.se`:

  standard error (continuous)

- `.r`:

  event count (discrete)

- `.n`:

  event count denominator (discrete, `agd_arm` only)

- `.E`:

  time at risk (discrete)

- `.Surv`:

  survival outcome of type
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) (time-to-event),
  nested by study arm

- `.sample_size`:

  sample size (`agd_*` only)

- `...`:

  other columns (typically covariates) from the original data frame

Objects of class `mlnmr_data` additionally have components:

- `n_int`:

  number of numerical integration points

- `int_names`:

  names of covariates with numerical integration points

- `int_cor`:

  correlation matrix for covariates used to generate numerical
  integration points

The `agd_arm` and `agd_contrast` tibbles have additional list columns
with prefix `.int_`, one for each covariate, which contain the numerical
integration points nested as length-`n_int` vectors within each row.

## See also

[`print.nma_data()`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_data.md)
for the print method displaying details of the network, and
[`plot.nma_data()`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_data.md)
for network plots.
