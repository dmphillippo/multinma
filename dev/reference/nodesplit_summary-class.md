# The `nodesplit_summary` class

The `nodesplit_summary` class contains posterior summary statistics for
node-splitting models, as a result of calling
[`summary()`](https://rdrr.io/r/base/summary.html) on a `nma_nodesplit`
or `nma_nodesplit_df` object.

## Details

Objects of class `nodesplit_summary` are tibble data frames, with one
row for each node-split comparison and columns:

- `trt1`, `trt2`:

  Treatments forming the comparison

- `summary`:

  A list column containing
  [nma_summary](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
  objects with the posterior summaries and draws for each of the
  node-splitting parameters

- `p_value`:

  Bayesian p-value for inconsistency

- `dic`:

  A list column containing
  [nma_dic](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-class.md)
  objects, giving the model fit statistics

The parameters included in `summary` are:

- `d_net`:

  Network estimate from the corresponding consistency model, if
  available

- `d_dir`:

  Direct estimate from the node-splitting model

- `d_ind`:

  Indirect estimate from the node-splitting model

- `omega`:

  Inconsistency factor \\\omega = d\_\mathrm{dir} - d\_\mathrm{ind}\\

- `tau`:

  Heterogeneity standard deviation from the node-splitting model, if a
  random effects model was fitted

- `tau_consistency`:

  Heterogeneity standard deviation from the corresponding consistency
  model, if available and if a random effects model was fitted
