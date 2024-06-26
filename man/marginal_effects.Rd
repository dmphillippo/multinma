% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_effects.R
\name{marginal_effects}
\alias{marginal_effects}
\title{Marginal treatment effects}
\usage{
marginal_effects(
  object,
  ...,
  mtype = c("difference", "ratio", "link"),
  all_contrasts = FALSE,
  trt_ref = NULL,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  predictive_distribution = FALSE,
  summary = TRUE
)
}
\arguments{
\item{object}{A \code{stan_nma} object created by \code{\link[=nma]{nma()}}.}

\item{...}{Arguments passed to \code{\link[=predict.stan_nma]{predict.stan_nma()}}, for example to specify
the covariate distribution and baseline risk for a target population, e.g.
\code{newdata}, \code{baseline}, and related arguments. For survival outcomes, \code{type}
can also be specified to determine the quantity from which to form a
marginal effect. For example, \code{type = "hazard"} with \code{mtype = "ratio"}
produces marginal hazard ratios, \code{type = "median"} with \code{mtype = "difference"} produces marginal median survival time differences, and so
on.}

\item{mtype}{The type of marginal effect to construct from the average
absolute effects, either \code{"difference"} (the default) for a difference of
absolute effects such as a risk difference, \code{"ratio"} for a ratio of
absolute effects such as a risk ratio, or \code{"link"} for a difference on the
scale of the link function used in fitting the model such as a marginal log
odds ratio.}

\item{all_contrasts}{Logical, generate estimates for all contrasts (\code{TRUE}),
or just the "basic" contrasts against the network reference treatment
(\code{FALSE})? Default \code{FALSE}.}

\item{trt_ref}{Reference treatment to construct relative effects against, if
\code{all_contrasts = FALSE}. By default, relative effects will be against the
network reference treatment. Coerced to character string.}

\item{probs}{Numeric vector of quantiles of interest to present in computed
summary, default \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}}

\item{predictive_distribution}{Logical, when a random effects model has been
fitted, should the predictive distribution for marginal effects in a new
study be returned? Default \code{FALSE}.}

\item{summary}{Logical, calculate posterior summaries? Default \code{TRUE}.}
}
\value{
A \link{nma_summary} object if \code{summary = TRUE}, otherwise a list
containing a 3D MCMC array of samples and (for regression models) a data
frame of study information.
}
\description{
Generate population-average marginal treatment effects. These are formed from
population-average absolute predictions, so this function is a wrapper around
\code{\link[=predict.stan_nma]{predict.stan_nma()}}.
}
\examples{
## Smoking cessation
\donttest{
# Run smoking RE NMA example if not already available
if (!exists("smk_fit_RE")) example("example_smk_re", run.donttest = TRUE)
}
\donttest{
# Marginal risk difference in each study population in the network
marginal_effects(smk_fit_RE, mtype = "difference")

# Since there are no covariates in the model, the marginal and conditional
# (log) odds ratios here coincide
marginal_effects(smk_fit_RE, mtype = "link")
relative_effects(smk_fit_RE)

# Marginal risk differences in a population with 67 observed events out of
# 566 individuals on No Intervention, corresponding to a Beta(67, 566 - 67)
# distribution on the baseline probability of response
(smk_rd_RE <- marginal_effects(smk_fit_RE,
                               baseline = distr(qbeta, 67, 566 - 67),
                               baseline_type = "response",
                               mtype = "difference"))
plot(smk_rd_RE)
}

## Plaque psoriasis ML-NMR
\donttest{
# Run plaque psoriasis ML-NMR example if not already available
if (!exists("pso_fit")) example("example_pso_mlnmr", run.donttest = TRUE)
}
\donttest{
# Population-average marginal probit differences in each study in the network
(pso_marg <- marginal_effects(pso_fit, mtype = "link"))
plot(pso_marg, ref_line = c(0, 1))

# Population-average marginal probit differences in a new target population,
# with means and SDs or proportions given by
new_agd_int <- data.frame(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1
)

# We need to add integration points to this data frame of new data
# We use the weighted mean correlation matrix computed from the IPD studies
new_agd_int <- add_integration(new_agd_int,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               cor = pso_net$int_cor,
                               n_int = 64)

# Population-average marginal probit differences of achieving PASI 75 in this
# target population, given a Normal(-1.75, 0.08^2) distribution on the
# baseline probit-probability of response on Placebo (at the reference levels
# of the covariates), are given by
(pso_marg_new <- marginal_effects(pso_fit,
                                  mtype = "link",
                                  newdata = new_agd_int,
                                  baseline = distr(qnorm, -1.75, 0.08)))
plot(pso_marg_new)
}

## Progression free survival with newly-diagnosed multiple myeloma
\donttest{
# Run newly-diagnosed multiple myeloma example if not already available
if (!exists("ndmm_fit")) example("example_ndmm", run.donttest = TRUE)
}
\donttest{
# We can produce a range of marginal effects from models with survival
# outcomes, specified with the mtype and type arguments. For example:

# Marginal survival probability difference at 5 years, all contrasts
marginal_effects(ndmm_fit, type = "survival", mtype = "difference",
                 times = 5, all_contrasts = TRUE)

# Marginal difference in RMST up to 5 years
marginal_effects(ndmm_fit, type = "rmst", mtype = "difference", times = 5)

# Marginal median survival time ratios
marginal_effects(ndmm_fit, type = "median", mtype = "ratio")

# Marginal log hazard ratios
# With no covariates in the model, these are constant over time and study
# populations, and are equal to the log hazard ratios from relative_effects()
plot(marginal_effects(ndmm_fit, type = "hazard", mtype = "link"),
     # The hazard is infinite at t=0 in some studies, giving undefined logHRs at t=0
     na.rm = TRUE)

# The NDMM vignette demonstrates the production of time-varying marginal
# hazard ratios from a ML-NMR model that includes covariates, see
# `vignette("example_ndmm")`

# Marginal survival difference over time
plot(marginal_effects(ndmm_fit, type = "survival", mtype = "difference"))
}

}
