# multinma 0.3.0.9000

* Feature: Node-splitting models for assessing inconsistency are now available
with `consistency = "nodesplit"` in `nma()`. Comparisons to split can be chosen
using the `nodesplit` argument, by default all possibly inconsistent comparisons
are chosen using `get_nodesplits()`. Node-splitting results can be summarised
with `summary.nma_nodesplit()` and plotted with `plot.nodesplit_summary()`.
* Feature: The correlation matrix for generating integration points with
`add_integration()` for ML-NMR models is now adjusted to the underlying Gaussian
copula, so that the output correlations of the integration points better match
the requested input correlations. A new argument `cor_adjust` controls this
behaviour, with options `"spearman"`, `"pearson"`, or `"none"`. Although these
correlations typically have little impact on the results, for strict
reproducibility the old behaviour from version 0.3.0 and below is available with
`cor_adjust = "legacy"`.
* Improvement: Factor order is now respected when `trt`, `study`, or `trt_class`
are factors, previously the order of levels was reset into natural sort order.
* Improvement: Update package website to Bootstrap 5 with release of pkgdown
2.0.0
* Fix: Model fitting is now robust to non-default settings of
`options("contrasts")`.
* Fix: `plot.nma_data()` no longer gives a ggplot deprecation warning (PR #6).
* Fix: Bug in `predict.stan_nma()` with a single covariate when `newdata` is a
`data.frame` (PR #7).
* Fix: Attempting to call `predict.stan_nma()` on a regression model with only
contrast data and no `newdata` or `baseline` specified now throws a descriptive
error message.

# multinma 0.3.0

* Feature: Added `baseline_type` and `baseline_level` arguments to
`predict.stan_nma()`, which allow baseline distributions to be specified on the
response or linear predictor scale, and at the individual or aggregate level.
* Feature: The `baseline` argument to `predict.stan_nma()` can now accept a
(named) list of baseline distributions if `newdata` contains multiple studies.
* Improvement: Misspecified `newdata` arguments to functions like
`relative_effects()` and `predict.stan_nma()` now give more informative error
messages.
* Fix: Constructing models with contrast-based data previously gave errors in
some scenarios (ML-NMR models, UME models, and in some cases AgD meta-regression
models).
* Fix: Ensure CRAN additional checks with `--run-donttest` run correctly.

# multinma 0.2.1

* Fix: Producing relative effect estimates for all contrasts using
`relative_effects()` with `all_contrasts = TRUE` no longer gives an error for
regression models.
* Fix: Specifying the covariate correlation matrix `cor` in `add_integration()`
is not required when only one covariate is present.
* Improvement: Added more detailed documentation on the likelihoods and link
functions available for each data type (`likelihood` and `link` arguments in
`nma()`).

# multinma 0.2.0

* Feature: The `set_*()` functions now accept `dplyr::mutate()` style semantics,
allowing inline variable transformations.
* Feature: Added ordered multinomial models, with helper function `multi()` for
specifying the outcomes. Accompanied by a new data set `hta_psoriasis` and
vignette.
* Feature: Implicit flat priors can now be specified, on any parameter, using
`flat()`.
* Improvement: `as.array.stan_nma()` is now much more efficient, meaning that
many post-estimation functions are also now much more efficient.
* Improvement: `plot.nma_dic()` is now more efficient, particularly with large
numbers of data points.
* Improvement: The layering of points when producing "dev-dev" plots using
`plot.nma_dic()` with multiple data types has been reversed for improved clarity
(now AgD over the top of IPD).
* Improvement: Aggregate-level predictions with `predict()` from ML-NMR / IPD
regression models are now calculated in a much more memory-efficient manner.
* Improvement: Added an overview of examples given in the vignettes.
* Improvement: Network plots with `weight_edges = TRUE` no longer produce
legends with non-integer values for the number of studies.
* Fix: `plot.nma_dic()` no longer gives an error when attempting to specify
`.width` argument when producing "dev-dev" plots.

# multinma 0.1.3

* Format DESCRIPTION to CRAN requirements

# multinma 0.1.2

* Wrapped long-running examples in `\donttest{}` instead of `\dontrun{}`

# multinma 0.1.1

* Reduced size of vignettes
* Added methods paper reference to DESCRIPTION
* Added zenodo DOI

# multinma 0.1.0

* Feature: Network plots, using a `plot()` method for `nma_data` objects.
* Feature: `as.igraph()`, `as_tbl_graph()` methods for `nma_data` objects.
* Feature: Produce relative effect estimates with `relative_effects()`,
posterior ranks with `posterior_ranks()`, and posterior rank probabilities with
`posterior_rank_probs()`. These will be study-specific when a regression model
is given.
* Feature: Produce predictions of absolute effects with a `predict()` method for
`stan_nma` objects.
* Feature: Plots of relative effects, ranks, predictions, and parameter
estimates via `plot.nma_summary()`.
* Feature: Optional `sample_size` argument for `set_agd_*()` that:
  - Enables centering of predictors (`center = TRUE`) in `nma()` when
    a regression model is given, replacing the `agd_sample_size` argument of `nma()`
  - Enables production of study-specific relative effects, rank probabilities,
    etc. for studies in the network when a regression model is given
  - Allows nodes in network plots to be weighted by sample size
* Feature: Plots of residual deviance contributions for a model and "dev-dev"
plots comparing residual deviance contributions between two models, using a
`plot()` method for `nma_dic` objects produced by `dic()`.
* Feature: Complementary log-log (cloglog) link function `link = "cloglog"` for
binomial likelihoods.
* Feature: Option to specify priors for heterogeneity on the standard deviation,
variance, or precision, with argument `prior_het_type`.
* Feature: Added log-Normal prior distribution.
* Feature: Plots of prior distributions vs. posterior distributions with
`plot_prior_posterior()`.
* Feature: Pairs plot method `pairs()`.
* Feature: Added vignettes with example analyses from the NICE TSDs and more.
* Fix: Random effects models with even moderate numbers of studies could be very
slow. These now run much more quickly, using a sparse representation of the RE
correlation matrix which is automatically enabled for sparsity above 90%
(roughly equivalent to 10 or more studies).

# multinma 0.0.1

* Initial release.
