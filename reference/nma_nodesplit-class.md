# The nma_nodesplit class

The `nma_nodesplit` and `nma_nodesplit_df` classes contains the results
from running a node-splitting model with the function
[`nma()`](https://dmphillippo.github.io/multinma/reference/nma.md).

## Details

Objects of class `nma_nodesplit` inherit from the
[stan_nma](https://dmphillippo.github.io/multinma/reference/stan_nma-class.md)
class, and contain the results of fitting a single node-split model.
They have one additional component, `nodesplit`, which gives the
comparison that was node-split as a length 2 vector.

Objects of class `nma_nodesplit_df` are tibble data frames with one row
for each node-split comparison and columns:

- `trt1`, `trt2`:

  Treatments forming the comparison

- `model`:

  A list column containing the results of each model as a
  `nma_nodesplit` object

Optionally, there will be an additional row for the consistency model if
this was fitted (e.g. by
`get_nodesplits(., include_consistency = TRUE)`) with `trt1` and `trt2`
both `NA`.
