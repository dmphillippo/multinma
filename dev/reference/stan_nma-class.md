# The stan_nma class

The `stan_nma` and `stan_mlnmr` classes contains the results from
running a model with the function
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md).

## Details

Objects of class `stan_nma` and `stan_mlnmr` have the following
components:

- `network`:

  The network data from which the model was run (class
  [nma_data](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  for `stan_nma`, or class
  [mlnmr_data](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  for `stan_mlnmr`)

- `stanfit`:

  The `stanfit` object returned by calling `sampling()` for the model

- `trt_effects`:

  Whether fixed or random effects were used (character string)

- `consistency`:

  The consistency/inconsistency model used (character string)

- `regression`:

  The regression model used (formula)

- `class_interactions`:

  If treatment classes and a regression model are specified, the model
  used for interactions within each class (common, exchangeable, or
  independent)

- `xbar`:

  A named vector of values used for centering

- `likelihood`:

  The likelihood used (character string)

- `link`:

  The link function used (character string)

- `priors`:

  A list containing the priors used (as
  [nma_prior](https://dmphillippo.github.io/multinma/dev/reference/nma_prior-class.md)
  objects)

- `basis`:

  For `mspline` and `pexp` models, a named list of spline bases for each
  study

The `stan_mlnmr` sub-class inherits from `stan_nma`, and differs only in
the class of the `network` object.
