# The nma_dic class

The `nma_dic` class contains details of the Deviance Information
Criterion (DIC), produced using the
[`dic()`](https://dmphillippo.github.io/multinma/reference/dic.md)
function.

## Details

Objects of class `nma_dic` have the following components:

- `dic`:

  The DIC value

- `pd`, `pv`:

  The effective number of parameters

- `resdev`:

  The total residual deviance

- `pointwise`:

  A list of data frames containing the pointwise contributions for the
  IPD and AgD.

- `resdev_array`:

  A 3D MCMC array \[Iterations, Chains, Parameters\] of posterior
  residual deviance samples.

## See also

[`dic()`](https://dmphillippo.github.io/multinma/reference/dic.md),
[`print.nma_dic()`](https://dmphillippo.github.io/multinma/reference/nma_dic-methods.md),
[`plot.nma_dic()`](https://dmphillippo.github.io/multinma/reference/plot.nma_dic.md).
