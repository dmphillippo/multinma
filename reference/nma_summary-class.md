# The `nma_summary` class

The `nma_summary` class contains posterior summary statistics of model
parameters or other quantities of interest, and the draws used to obtain
these statistics.

## Details

Objects of class `nma_summary` have the following components:

- summary:

  A data frame containing the computed summary statistics. Column `.trt`
  indicates the corresponding treatment, or columns `.trta` and `.trtb`
  indicate the corresponding contrast (`.trtb` vs. `.trta`). If a
  regression model was fitted with effect modifier interactions with
  treatment, these summaries will be study-specific. In this case, the
  corresponding study population is indicated in the `.study` column. If
  a multinomial model was fitted, the `.category` column indicates the
  corresponding category.

- sims:

  A 3D array \[Iteration, Chain, Parameter\] of MCMC simulations

- studies:

  (Optional) A data frame containing study information, printed along
  with the corresponding summary statistics if `summary` contains a
  `.study` column. Should have a matching `.study` column.

The following attributes may also be set:

- xlab:

  Label for x axis in plots, usually either `"Treatment"` or
  `"Contrast"`.

- ylab:

  Label for y axis in plots, usually used for the scale e.g.
  `"log Odds Ratio"`.

The subclass `nma_rank_probs` is used by the function
[`posterior_rank_probs()`](https://dmphillippo.github.io/multinma/reference/posterior_ranks.md),
and contains posterior rank probabilities. This subclass does not have a
`sims` component, as the rank probabilities are themselves posterior
summaries of the ranks (i.e. they do not have a posterior distribution).
The posterior ranks from which the rank probabilities are calculated may
be obtained from
[`posterior_ranks()`](https://dmphillippo.github.io/multinma/reference/posterior_ranks.md).
