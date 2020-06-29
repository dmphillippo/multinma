# multinma 0.1.3.9000

# multinma 0.1.3

* Format DESCRIPTION to CRAN requirements

# multinma 0.1.2

* Wrapped long-running examples in \donttest{} instead of \dontrun{}

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
