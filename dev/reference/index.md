# Package index

## Package overview

- [`multinma-package`](https://dmphillippo.github.io/multinma/dev/reference/multinma-package.md)
  [`multinma`](https://dmphillippo.github.io/multinma/dev/reference/multinma-package.md)
  : multinma: A Package for Network Meta-Analysis of Individual and
  Aggregate Data in Stan

## Defining a network

Setting up a network from different data sources, creating network
plots.

- [`set_agd_arm()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_arm.md)
  : Set up arm-based aggregate data

- [`set_agd_contrast()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_contrast.md)
  : Set up contrast-based aggregate data

- [`set_agd_surv()`](https://dmphillippo.github.io/multinma/dev/reference/set_agd_surv.md)
  : Set up aggregate survival data

- [`set_ipd()`](https://dmphillippo.github.io/multinma/dev/reference/set_ipd.md)
  : Set up individual patient data

- [`combine_network()`](https://dmphillippo.github.io/multinma/dev/reference/combine_network.md)
  : Combine multiple data sources into one network

- [`multi()`](https://dmphillippo.github.io/multinma/dev/reference/multi.md)
  : Multinomial outcome data

- [`print(`*`<nma_data>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_data.md)
  [`print(`*`<mlnmr_data>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_data.md)
  :

  Print `nma_data` objects

- [`plot(`*`<nma_data>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_data.md)
  : Network plots

- [`as.igraph(`*`<nma_data>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/graph_conversion.md)
  [`as_tbl_graph(`*`<nma_data>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/graph_conversion.md)
  : Convert networks to graph objects

- [`nma_data-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  [`nma_data`](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  [`mlnmr_data`](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  [`mlnmr_data-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_data-class.md)
  : The nma_data class

- [`is_network_connected()`](https://dmphillippo.github.io/multinma/dev/reference/is_network_connected.md)
  : Check network connectedness

## Setting up numerical integration (ML-NMR only)

Multilevel network meta-regression models require numerical integration
points to be specified for the distributions of covariates in each
aggregate data study in the network.

- [`add_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md)
  [`unnest_integration()`](https://dmphillippo.github.io/multinma/dev/reference/add_integration.md)
  : Add numerical integration points to aggregate data
- [`distr()`](https://dmphillippo.github.io/multinma/dev/reference/distr.md)
  : Specify a general marginal distribution
- [`qbern()`](https://dmphillippo.github.io/multinma/dev/reference/Bernoulli.md)
  [`pbern()`](https://dmphillippo.github.io/multinma/dev/reference/Bernoulli.md)
  [`dbern()`](https://dmphillippo.github.io/multinma/dev/reference/Bernoulli.md)
  : The Bernoulli Distribution
- [`qgamma()`](https://dmphillippo.github.io/multinma/dev/reference/GammaDist.md)
  [`dgamma()`](https://dmphillippo.github.io/multinma/dev/reference/GammaDist.md)
  [`pgamma()`](https://dmphillippo.github.io/multinma/dev/reference/GammaDist.md)
  : The Gamma distribution
- [`dgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  [`pgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  [`qgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  : Generalised Student's t distribution (with location and scale)
- [`dlogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  [`plogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  [`qlogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  : Log Student's t distribution
- [`dlogitnorm()`](https://dmphillippo.github.io/multinma/dev/reference/logitNormal.md)
  [`plogitnorm()`](https://dmphillippo.github.io/multinma/dev/reference/logitNormal.md)
  [`qlogitnorm()`](https://dmphillippo.github.io/multinma/dev/reference/logitNormal.md)
  : The logit Normal distribution

## Prior distributions

Specify and summarise prior distributions.

- [`normal()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`half_normal()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`log_normal()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`cauchy()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`half_cauchy()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`student_t()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`half_student_t()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`log_student_t()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`exponential()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  [`flat()`](https://dmphillippo.github.io/multinma/dev/reference/priors.md)
  : Prior distributions
- [`summary(`*`<nma_prior>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.nma_prior.md)
  : Summary of prior distributions
- [`nma_prior-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_prior-class.md)
  [`nma_prior`](https://dmphillippo.github.io/multinma/dev/reference/nma_prior-class.md)
  : The nma_prior class
- [`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/dev/reference/plot_prior_posterior.md)
  : Plot prior vs posterior distribution
- [`dgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  [`pgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  [`qgent()`](https://dmphillippo.github.io/multinma/dev/reference/generalised_t.md)
  : Generalised Student's t distribution (with location and scale)
- [`dlogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  [`plogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  [`qlogt()`](https://dmphillippo.github.io/multinma/dev/reference/log_t.md)
  : Log Student's t distribution

## Model fitting

Model specification and fitting is accomplished using the
[`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
function.

- [`nma()`](https://dmphillippo.github.io/multinma/dev/reference/nma.md)
  : Network meta-analysis models

- [`print(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/print.stan_nma.md)
  :

  Print `stan_nma` objects

- [`summary(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.stan_nma.md)
  [`plot(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.stan_nma.md)
  :

  Posterior summaries from `stan_nma` objects

- [`pairs(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/pairs.stan_nma.md)
  :

  Matrix of plots for a `stan_nma` object

- [`stan_nma-class`](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
  [`stan_nma`](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
  [`stan_mlnmr`](https://dmphillippo.github.io/multinma/dev/reference/stan_nma-class.md)
  : The stan_nma class

- [`adapt_delta`](https://dmphillippo.github.io/multinma/dev/reference/adapt_delta.md)
  : Target average acceptance probability

- [`RE_cor()`](https://dmphillippo.github.io/multinma/dev/reference/random_effects.md)
  [`which_RE()`](https://dmphillippo.github.io/multinma/dev/reference/random_effects.md)
  : Random effects structure

- [`.default()`](https://dmphillippo.github.io/multinma/dev/reference/default_values.md)
  [`.is_default()`](https://dmphillippo.github.io/multinma/dev/reference/default_values.md)
  : Set default values

## Model checking and comparison

Checking model fit and comparing models.

- [`plot_prior_posterior()`](https://dmphillippo.github.io/multinma/dev/reference/plot_prior_posterior.md)
  : Plot prior vs posterior distribution

- [`plot_integration_error()`](https://dmphillippo.github.io/multinma/dev/reference/plot_integration_error.md)
  : Plot numerical integration error

- [`dic()`](https://dmphillippo.github.io/multinma/dev/reference/dic.md)
  : Deviance Information Criterion (DIC)

- [`print(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  [`as.data.frame(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  [`as.tibble(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  [`as_tibble(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  [`as.array(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  [`as.matrix(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-methods.md)
  :

  Methods for `nma_dic` objects

- [`plot(`*`<nma_dic>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_dic.md)
  : Plots of model fit diagnostics

- [`nma_dic-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-class.md)
  [`nma_dic`](https://dmphillippo.github.io/multinma/dev/reference/nma_dic-class.md)
  : The nma_dic class

- [`loo(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/loo.md)
  [`waic(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/loo.md)
  :

  Model comparison using the `loo` package

## Node-splitting

Generate and summarise node-splitting models for assessing
inconsistency.

- [`get_nodesplits()`](https://dmphillippo.github.io/multinma/dev/reference/get_nodesplits.md)
  [`has_direct()`](https://dmphillippo.github.io/multinma/dev/reference/get_nodesplits.md)
  [`has_indirect()`](https://dmphillippo.github.io/multinma/dev/reference/get_nodesplits.md)
  : Direct and indirect evidence

- [`nma_nodesplit-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_nodesplit-class.md)
  [`nma_nodesplit`](https://dmphillippo.github.io/multinma/dev/reference/nma_nodesplit-class.md)
  [`nma_nodesplit_df`](https://dmphillippo.github.io/multinma/dev/reference/nma_nodesplit-class.md)
  [`nma_nodesplit_df-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_nodesplit-class.md)
  : The nma_nodesplit class

- [`print(`*`<nma_nodesplit_df>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_nodesplit_df.md)
  [`print(`*`<nma_nodesplit>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/print.nma_nodesplit_df.md)
  :

  Print `nma_nodesplit_df` objects

- [`summary(`*`<nma_nodesplit_df>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.nma_nodesplit_df.md)
  [`summary(`*`<nma_nodesplit>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.nma_nodesplit_df.md)
  [`plot(`*`<nma_nodesplit>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.nma_nodesplit_df.md)
  [`plot(`*`<nma_nodesplit_df>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/summary.nma_nodesplit_df.md)
  : Summarise the results of node-splitting models

- [`nodesplit_summary-class`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-class.md)
  [`nodesplit_summary`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-class.md)
  :

  The `nodesplit_summary` class

- [`print(`*`<nodesplit_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-methods.md)
  [`as_tibble(`*`<nodesplit_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-methods.md)
  [`as.tibble.nodesplit_summary()`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-methods.md)
  [`as.data.frame(`*`<nodesplit_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nodesplit_summary-methods.md)
  :

  Methods for `nodesplit_summary` objects

- [`plot(`*`<nodesplit_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nodesplit_summary.md)
  : Plots of node-splitting models

## Posterior summaries and working with fitted models

Producing and plotting relative effects, absolute predictions, marginal
effects, posterior ranks and rank probabilities. Converting to MCMC
arrays and matrices.

- [`relative_effects()`](https://dmphillippo.github.io/multinma/dev/reference/relative_effects.md)
  : Relative treatment effects

- [`marginal_effects()`](https://dmphillippo.github.io/multinma/dev/reference/marginal_effects.md)
  : Marginal treatment effects

- [`predict(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md)
  [`predict(`*`<stan_nma_surv>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/predict.stan_nma.md)
  : Predictions of absolute effects from NMA models

- [`posterior_ranks()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md)
  [`posterior_rank_probs()`](https://dmphillippo.github.io/multinma/dev/reference/posterior_ranks.md)
  : Treatment rankings and rank probabilities

- [`print(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.data.frame(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.tibble(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as_tibble(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.array(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.matrix(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.array(`*`<nma_rank_probs>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  [`as.matrix(`*`<nma_rank_probs>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-methods.md)
  :

  Methods for `nma_summary` objects

- [`plot(`*`<nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md)
  [`plot(`*`<nma_parameter_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md)
  [`plot(`*`<nma_rank_probs>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md)
  [`plot(`*`<surv_nma_summary>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/plot.nma_summary.md)
  : Plots of summary results

- [`nma_summary-class`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
  [`nma_summary`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
  [`nma_rank_probs`](https://dmphillippo.github.io/multinma/dev/reference/nma_summary-class.md)
  :

  The `nma_summary` class

- [`as.array(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/as.array.stan_nma.md)
  [`as.data.frame(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/as.array.stan_nma.md)
  [`as_tibble(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/as.array.stan_nma.md)
  [`as.tibble(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/as.array.stan_nma.md)
  [`as.matrix(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/as.array.stan_nma.md)
  : Convert samples into arrays, matrices, or data frames

- [`summary(`*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  [`print(`*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  [`plot(`*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  [`names(`*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  [`` `names<-`( ``*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/mcmc_array-class.md)
  : Working with 3D MCMC arrays

- [`as.stanfit()`](https://dmphillippo.github.io/multinma/dev/reference/as.stanfit.md)
  : as.stanfit

- [`bind_chains()`](https://dmphillippo.github.io/multinma/dev/reference/bind_chains.md)
  [`cbind(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/bind_chains.md)
  [`cbind(`*`<mcmc_array>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/bind_chains.md)
  : Bind chains

## M-spline hazards

Functions for flexibile M-splines on the baseline hazard.

- [`dmspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  [`pmspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  [`qmspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  [`hmspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  [`Hmspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  [`rmst_mspline()`](https://dmphillippo.github.io/multinma/dev/reference/mspline.md)
  : Distribution functions for M-spline baseline hazards
- [`make_knots()`](https://dmphillippo.github.io/multinma/dev/reference/make_knots.md)
  : Knot locations for M-spline baseline hazard models
- [`knots(`*`<stan_nma>`*`)`](https://dmphillippo.github.io/multinma/dev/reference/knots.stan_nma.md)
  : Knot locations for a fitted model
- [`softmax()`](https://dmphillippo.github.io/multinma/dev/reference/softmax.md)
  [`inv_softmax()`](https://dmphillippo.github.io/multinma/dev/reference/softmax.md)
  : Softmax transform

## ggplot functions

Functions for creating or customising ggplot outputs.

- [`theme_multinma()`](https://dmphillippo.github.io/multinma/dev/reference/theme_multinma.md)
  : Plot theme for multinma plots
- [`geom_km()`](https://dmphillippo.github.io/multinma/dev/reference/geom_km.md)
  : Kaplan-Meier curves of survival data

## Datasets

Datasets used for examples and vignettes.

- [`atrial_fibrillation`](https://dmphillippo.github.io/multinma/dev/reference/atrial_fibrillation.md)
  : Stroke prevention in atrial fibrillation patients
- [`bcg_vaccine`](https://dmphillippo.github.io/multinma/dev/reference/bcg_vaccine.md)
  : BCG vaccination
- [`blocker`](https://dmphillippo.github.io/multinma/dev/reference/blocker.md)
  : Beta blockers to prevent mortality after MI
- [`certolizumab`](https://dmphillippo.github.io/multinma/dev/reference/certolizumab.md)
  : Certolizumab
- [`diabetes`](https://dmphillippo.github.io/multinma/dev/reference/diabetes.md)
  : Incidence of diabetes in trials of antihypertensive drugs
- [`dietary_fat`](https://dmphillippo.github.io/multinma/dev/reference/dietary_fat.md)
  : Reduced dietary fat to prevent mortality
- [`hta_psoriasis`](https://dmphillippo.github.io/multinma/dev/reference/hta_psoriasis.md)
  : HTA Plaque Psoriasis
- [`ndmm_ipd`](https://dmphillippo.github.io/multinma/dev/reference/ndmm.md)
  [`ndmm_agd`](https://dmphillippo.github.io/multinma/dev/reference/ndmm.md)
  [`ndmm_agd_covs`](https://dmphillippo.github.io/multinma/dev/reference/ndmm.md)
  : Newly diagnosed multiple myeloma
- [`parkinsons`](https://dmphillippo.github.io/multinma/dev/reference/parkinsons.md)
  : Mean off-time reduction in Parkison's disease
- [`plaque_psoriasis_ipd`](https://dmphillippo.github.io/multinma/dev/reference/plaque_psoriasis.md)
  [`plaque_psoriasis_agd`](https://dmphillippo.github.io/multinma/dev/reference/plaque_psoriasis.md)
  : Plaque psoriasis data
- [`smoking`](https://dmphillippo.github.io/multinma/dev/reference/smoking.md)
  : Smoking cessation data
- [`social_anxiety`](https://dmphillippo.github.io/multinma/dev/reference/social_anxiety.md)
  : Social Anxiety
- [`statins`](https://dmphillippo.github.io/multinma/dev/reference/statins.md)
  : Statins for cholesterol lowering
- [`thrombolytics`](https://dmphillippo.github.io/multinma/dev/reference/thrombolytics.md)
  : Thrombolytic treatments data
- [`transfusion`](https://dmphillippo.github.io/multinma/dev/reference/transfusion.md)
  : Granulocyte transfusion in patients with neutropenia or neutrophil
  dysfunction
- [`example_ndmm`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_ndmm.md)
  : Example newly-diagnosed multiple myeloma
- [`example_pso_mlnmr`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_pso_mlnmr.md)
  : Example plaque psoriasis ML-NMR
- [`example_smk_fe`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_smk_fe.md)
  : Example smoking FE NMA
- [`example_smk_nodesplit`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_smk_nodesplit.md)
  : Example smoking node-splitting
- [`example_smk_re`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_smk_re.md)
  : Example smoking RE NMA
- [`example_smk_ume`](https://dmphillippo.github.io/multinma/dev/reference/aa_example_smk_ume.md)
  : Example smoking UME NMA
