% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nma.R
\name{nma}
\alias{nma}
\title{Network meta-analysis models}
\usage{
nma(
  network,
  consistency = c("consistency", "ume", "nodesplit"),
  trt_effects = c("fixed", "random"),
  regression = NULL,
  class_interactions = c("common", "exchangeable", "independent"),
  class_effects = c("independent", "common", "exchangeable"),
  class_sd = c("independent", "common"),
  likelihood = NULL,
  link = NULL,
  ...,
  nodesplit = get_nodesplits(network, include_consistency = TRUE),
  prior_intercept = .default(normal(scale = 100)),
  prior_trt = .default(normal(scale = 10)),
  prior_het = .default(half_normal(scale = 5)),
  prior_het_type = c("sd", "var", "prec"),
  prior_reg = .default(normal(scale = 10)),
  prior_aux = .default(),
  prior_aux_reg = .default(),
  prior_class_mean = .default(normal(scale = 10)),
  prior_class_sd = .default(half_normal(scale = 5)),
  aux_by = NULL,
  aux_regression = NULL,
  QR = FALSE,
  center = TRUE,
  adapt_delta = NULL,
  int_thin = 0,
  int_check = TRUE,
  mspline_degree = 3,
  n_knots = 7,
  knots = NULL,
  mspline_basis = NULL
)
}
\arguments{
\item{network}{An \code{nma_data} object, as created by the functions \verb{set_*()},
\code{combine_network()}, or \code{add_integration()}}

\item{consistency}{Character string specifying the type of (in)consistency
model to fit, either \code{"consistency"}, \code{"ume"}, or \code{"nodesplit"}}

\item{trt_effects}{Character string specifying either \code{"fixed"} or \code{"random"}
effects}

\item{regression}{A one-sided model formula, specifying the prognostic and
effect-modifying terms for a regression model. Any references to treatment
should use the \code{.trt} special variable, for example specifying effect
modifier interactions as \code{variable:.trt} (see details).}

\item{class_interactions}{Character string specifying whether effect modifier
interactions are specified as \code{"common"}, \code{"exchangeable"}, or
\code{"independent"}.}

\item{class_effects}{Character string specifying a model for treatment class effects,
either \code{"independent"} (the default), \code{"exchangeable"}, or \code{"common"}.}

\item{class_sd}{Character string specifying whether the class standard deviations in a
class effects model should be \code{"independent"} (i.e. separate for each class, the default),
or \code{"common"} (i.e. shared across all classes). Alternatively this can be a list of
character vectors, each of which describe a set classes for which to share a common class SD;
any list names will be used to name the output parameters, otherwise the name will be taken
from the first class in each set.}

\item{likelihood}{Character string specifying a likelihood, if unspecified
will be inferred from the data (see details)}

\item{link}{Character string specifying a link function, if unspecified will
default to the canonical link (see details)}

\item{...}{Further arguments passed to
\code{\link[rstan:stanmodel-method-sampling]{sampling()}}, such as \code{iter},
\code{chains}, \code{cores}, etc.}

\item{nodesplit}{For \code{consistency = "nodesplit"}, the comparison(s) to split
in the node-splitting model(s). Either a length 2 vector giving the
treatments in a single comparison, or a 2 column data frame listing
multiple treatment comparisons to split in turn. By default, all possible
comparisons will be chosen (see \code{\link[=get_nodesplits]{get_nodesplits()}}).}

\item{prior_intercept}{Specification of prior distribution for the intercept}

\item{prior_trt}{Specification of prior distribution for the treatment
effects}

\item{prior_het}{Specification of prior distribution for the heterogeneity
(if \code{trt_effects = "random"})}

\item{prior_het_type}{Character string specifying whether the prior
distribution \code{prior_het} is placed on the heterogeneity standard deviation
\eqn{\tau} (\code{"sd"}, the default), variance \eqn{\tau^2} (\code{"var"}), or
precision \eqn{1/\tau^2} (\code{"prec"}).}

\item{prior_reg}{Specification of prior distribution for the regression
coefficients (if \code{regression} formula specified)}

\item{prior_aux}{Specification of prior distribution for the auxiliary
parameter, if applicable (see details). For \code{likelihood = "gengamma"} this
should be a list of prior distributions with elements \code{sigma} and \code{k}.}

\item{prior_aux_reg}{Specification of prior distribution for the auxiliary
regression parameters, if \code{aux_regression} is specified (see details).}

\item{prior_class_mean}{Specification of prior distribution for the
treatment class means (if \code{class_effects = "exchangeable"}).}

\item{prior_class_sd}{Specification of prior distribution for the
treatment class standard deviations (if \code{class_effects = "exchangeable"}).}

\item{aux_by}{Vector of variable names listing the variables to stratify the
auxiliary parameters by. Currently only used for survival models, see
details. Cannot be used with \code{aux_regression}.}

\item{aux_regression}{A one-sided model formula giving a regression model for
the auxiliary parameters. Currently only used for survival models, see
details. Cannot be used with \code{aux_by}.}

\item{QR}{Logical scalar (default \code{FALSE}), whether to apply a QR
decomposition to the model design matrix}

\item{center}{Logical scalar (default \code{TRUE}), whether to center the
(numeric) regression terms about the overall means}

\item{adapt_delta}{See \link{adapt_delta} for details}

\item{int_thin}{A single integer value, the thinning factor for returning
cumulative estimates of integration error. Saving cumulative estimates is
disabled by \code{int_thin = 0}, which is the default.}

\item{int_check}{Logical, check sufficient accuracy of numerical integration
by fitting half of the chains with \code{n_int/2}? When \code{TRUE}, \code{Rhat} and
\code{n_eff} diagnostic warnings will be given if numerical integration has not
sufficiently converged (suggesting increasing \code{n_int} in
\code{\link[=add_integration]{add_integration()}}). Default \code{TRUE}, but disabled (\code{FALSE}) when
\code{int_thin > 0}.}

\item{mspline_degree}{Non-negative integer giving the degree of the M-spline
polynomial for \code{likelihood = "mspline"}. Piecewise exponential hazards
(\code{likelihood = "pexp"}) are a special case with \code{mspline_degree = 0}.}

\item{n_knots}{For \code{mspline} and \code{pexp} likelihoods, a non-negative integer
giving the number of internal knots for partitioning the baseline hazard
into intervals. The knot locations within each study will be determined by
the corresponding quantiles of the observed event times, plus boundary
knots at the earliest entry time (0 with no delayed entry) and the maximum
event/censoring time. For example, with \code{n_knots = 3}, the internal knot
locations will be at the 25\%, 50\%, and 75\% quantiles of the observed event
times. The default is \code{n_knots = 7}; overfitting is avoided by shrinking
towards a constant hazard with a random walk prior (see details). If
\code{aux_regression} is specified then a single set of knot locations will be
calculated across all studies in the network. See \code{\link[=make_knots]{make_knots()}} for more
details on the knot positioning algorithms. Ignored when \code{knots} is
specified.}

\item{knots}{For \code{mspline} and \code{pexp} likelihoods, a named list of numeric
vectors of knot locations (including boundary knots) for each of the
studies in the network. Currently, each vector must have the same length
(i.e. each study must use the same number of knots). Alternatively, a
single numeric vector of knot locations can be provided which will be
shared across all studies in the network. If unspecified (the default), the
knots will be chosen based on \code{n_knots} as described above. If
\code{aux_regression} is specified then \code{knots} should be a single numeric
vector of knot locations which will be shared across all studies in the
network. \code{\link[=make_knots]{make_knots()}} can be used to help specify \code{knots} directly, or to
investigate knot placement prior to model fitting.}

\item{mspline_basis}{Instead of specifying \code{mspline_degree} and \code{n_knots} or
\code{knots}, a named list of M-spline bases (one for each study) can be
provided with \code{mspline_basis} which will be used directly. In this case,
all other M-spline options will be ignored.}
}
\value{
\code{nma()} returns a \link{stan_nma} object, except when \code{consistency = "nodesplit"} when a \link{nma_nodesplit} or \link{nma_nodesplit_df} object is
returned. \code{nma.fit()} returns a \code{\link[rstan:stanfit-class]{stanfit}}
object.
}
\description{
The \code{nma} function fits network meta-analysis and (multilevel) network
meta-regression models in Stan.
}
\details{
When specifying a model formula in the \code{regression} argument, the
usual formula syntax is available (as interpreted by \code{\link[=model.matrix]{model.matrix()}}). The
only additional requirement here is that the special variable \code{.trt} should
be used to refer to treatment. For example, effect modifier interactions
should be specified as \code{variable:.trt}. Prognostic (main) effects and
interactions can be included together compactly as \code{variable*.trt}, which
expands to \code{variable + variable:.trt} (plus \code{.trt}, which is already in the
NMA model).

For the advanced user, the additional specials \code{.study} and \code{.trtclass} are
also available, and refer to studies and (if specified) treatment classes
respectively. When node-splitting models are fitted (\code{consistency = "nodesplit"}) the special \code{.omega} is available, indicating the arms to
which the node-splitting inconsistency factor is added.

See \code{\link[multinma:priors]{?priors}} for details on prior
specification. Default prior distributions are available, but may not be
appropriate for the particular setting and will raise a warning if used. No
attempt is made to tailor these defaults to the data provided. Please
consider appropriate prior distributions for the particular setting,
accounting for the scales of outcomes and covariates, etc. The function
\code{\link[=plot_prior_posterior]{plot_prior_posterior()}} may be useful in examining the influence of the
chosen prior distributions on the posterior distributions, and the
\code{\link[multinma:summary.nma_prior]{summary()}} method for \code{nma_prior}
objects prints prior intervals.
}
\section{Likelihoods and link functions}{

Currently, the following likelihoods and link functions are supported for
each data type:\tabular{lll}{
   \strong{Data type} \tab \strong{Likelihood} \tab \strong{Link function} \cr
   \strong{Binary} \tab \code{bernoulli}, \code{bernoulli2} \tab \code{logit}, \code{probit}, \code{cloglog} \cr
   \strong{Count} \tab \code{binomial}, \code{binomial2} \tab \code{logit}, \code{probit}, \code{cloglog} \cr
   \strong{Rate} \tab \code{poisson} \tab \code{log} \cr
   \strong{Continuous} \tab \code{normal} \tab \code{identity}, \code{log} \cr
   \strong{Ordered} \tab \code{ordered} \tab \code{logit}, \code{probit}, \code{cloglog} \cr
   \strong{Survival} \tab \code{exponential}, \code{weibull}, \code{gompertz}, \code{exponential-aft}, \code{weibull-aft}, \code{lognormal}, \code{loglogistic}, \code{gamma}, \code{gengamma}, \code{mspline}, \code{pexp} \tab \code{log} \cr
}


The \code{bernoulli2} and \code{binomial2} likelihoods correspond to a two-parameter
Binomial likelihood for arm-based AgD, which more closely matches the
underlying Poisson Binomial distribution for the summarised aggregate
outcomes in a ML-NMR model than the typical (one parameter) Binomial
distribution \insertCite{@see @methods_paper}{multinma}.

When a \code{cloglog} link is used, including an offset for log follow-up time
(i.e. \code{regression = ~offset(log(time))}) results in a model on the log
hazard \insertCite{@see @TSD2}{multinma}.

For survival data, all accelerated failure time models (\code{exponential-aft},
\code{weibull-aft}, \code{lognormal}, \code{loglogistic}, \code{gamma}, \code{gengamma}) are
parameterised so that the treatment effects and any regression parameters
are log Survival Time Ratios (i.e. a coefficient of \eqn{\log(2)} means
that the treatment or covariate is associated with a doubling of expected
survival time). These can be converted to log Acceleration Factors using
the relation \eqn{\log(\mathrm{AF}) = -\log(\mathrm{STR})} (or equivalently
\eqn{\mathrm{AF} = 1/\mathrm{STR}}).

Further details on each likelihood and link function are given by
\insertCite{TSD2;textual}{multinma}.
}

\section{Auxiliary parameters}{

Auxiliary parameters are only present in the following models.
\subsection{Normal likelihood with IPD}{

When a Normal likelihood is fitted to IPD, the auxiliary parameters are the
arm-level standard deviations \eqn{\sigma_{jk}} on treatment \eqn{k} in
study \eqn{j}.
}

\subsection{Ordered multinomial likelihood}{

When fitting a model to \eqn{M} ordered outcomes, the auxiliary parameters
are the latent cutoffs between each category, \eqn{c_0 < c_1 < \dots <
  c_M}. Only \eqn{c_2} to \eqn{c_{M-1}} are estimated; we fix \eqn{c_0 =
  -\infty}, \eqn{c_1 = 0}, and \eqn{c_M = \infty}. When specifying priors for
these latent cutoffs, we choose to specify priors on the \emph{differences}
\eqn{c_{m+1} - c_m}. Stan automatically truncates any priors so that the
ordering constraints are satisfied.
}

\subsection{Survival (time-to-event) likelihoods}{

All survival likelihoods except the \code{exponential} and \code{exponential-aft}
likelihoods have auxiliary parameters. These are typically study-specific
shape parameters \eqn{\gamma_j>0}, except for the \code{lognormal} likelihood
where the auxiliary parameters are study-specific standard deviations on
the log scale \eqn{\sigma_j>0}.

The \code{gengamma} likelihood has two sets of auxiliary parameters,
study-specific scale parameters \eqn{\sigma_j>0} and shape parameters
\eqn{k_j}, following the parameterisation of
\insertCite{Lawless1980;textual}{multinma}, which permits a range of
behaviours for the baseline hazard including increasing, decreasing,
bathtub and arc-shaped hazards. This parameterisation is related to that
discussed by \insertCite{Cox2007;textual}{multinma} and implemented in the
\code{flexsurv} package with \eqn{Q = k^{-0.5}}. The parameterisation used here
effectively bounds the shape parameter \eqn{k} away from numerical
instabilities as \eqn{k \rightarrow \infty} (i.e. away from \eqn{Q
  \rightarrow 0}, the log-Normal distribution) via the prior distribution.
Implicitly, this parameterisation is restricted to \eqn{Q > 0} and so
certain survival distributions like the inverse-Gamma and inverse-Weibull
are not part of the parameter space; however, \eqn{Q > 0} still encompasses
the other survival distributions implemented in this package.

For the \code{mspline} and \code{pexp} likelihoods, the auxiliary parameters are the
spline coefficients for each study. These form a unit simplex (i.e. lie
between 0 and 1, and sum to 1), and are given a random walk prior
distribution. \code{prior_aux} specifies the hyperprior on the random walk
standard deviation \eqn{\sigma} which controls the level of smoothing of
the baseline hazard, with \eqn{\sigma = 0} corresponding to a constant
baseline hazard.

The auxiliary parameters can be stratified by additional factors through
the \code{aux_by} argument. For example, to allow the shape of the baseline
hazard to vary between treatment arms as well as studies, use \code{aux_by = c(".study", ".trt")}. (Technically, \code{.study} is always included in the
stratification even if omitted from \code{aux_by}, but we choose here to make
the stratification explicit.) This is a common way of relaxing the
proportional hazards assumption. The default is equivalent to \code{aux_by = ".study"} which stratifies the auxiliary parameters by study, as described
above.

A regression model may be specified on the auxiliary parameters using
\code{aux_regression}. This is useful if we wish to model departures from
non-proportionality, rather than allowing the baseline hazards to be
completely independent using \code{aux_by}. This is necessary if absolute
predictions (e.g. survival curves) are required in a population for
unobserved combinations of covariates; for example, if \code{aux_by = .trt} then
absolute predictions may only be produced for the observed treatment arms
in each study population, whereas if \code{aux_regression = ~.trt} then absolute
predictions can be produced for all treatments in any population. For
\code{mspline} and \code{pexp} likelihoods, the regression coefficients are smoothed
over time using a random walk prior to avoid overfitting: \code{prior_aux_reg}
specifies the hyperprior for the random walk standard deviation. For other
parametric likelihoods, \code{prior_aux_reg} specifies the prior for the
auxiliary regression coefficients.
}
}

\examples{
## Smoking cessation NMA
# Set up network of smoking cessation data
head(smoking)

smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Print details
smk_net

\donttest{
# Fitting a fixed effect model
smk_fit_FE <- nma(smk_net,
                  trt_effects = "fixed",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100))

smk_fit_FE
}

\donttest{
# Fitting a random effects model
smk_fit_RE <- nma(smk_net, 
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = normal(scale = 5))

smk_fit_RE
}

\donttest{
# Fitting an unrelated mean effects (inconsistency) model
smk_fit_RE_UME <- nma(smk_net, 
                      consistency = "ume",
                      trt_effects = "random",
                      prior_intercept = normal(scale = 100),
                      prior_trt = normal(scale = 100),
                      prior_het = normal(scale = 5))

smk_fit_RE_UME
}

\donttest{
# Fitting all possible node-splitting models
smk_fit_RE_nodesplit <- nma(smk_net, 
                            consistency = "nodesplit",
                            trt_effects = "random",
                            prior_intercept = normal(scale = 100),
                            prior_trt = normal(scale = 100),
                            prior_het = normal(scale = 5))
}

\donttest{
# Summarise the node-splitting results
summary(smk_fit_RE_nodesplit)
}

## Plaque psoriasis ML-NMR
# Set up plaque psoriasis network combining IPD and AgD
library(dplyr)
pso_ipd <- filter(plaque_psoriasis_ipd,
                  studyc \%in\% c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"))

pso_agd <- filter(plaque_psoriasis_agd,
                  studyc == "FIXTURE")

head(pso_ipd)
head(pso_agd)

pso_ipd <- pso_ipd \%>\%
  mutate(# Variable transformations
    bsa = bsa / 100,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    weight = weight / 10,
    durnpso = durnpso / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn \%in\% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker"),
    # Check complete cases for covariates of interest
    complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  )

pso_agd <- pso_agd \%>\%
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    prevsys = prevsys / 100,
    psa = psa / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn \%in\% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker")
  )

# Exclude small number of individuals with missing covariates
pso_ipd <- filter(pso_ipd, complete)

pso_net <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = pasi75,
          trt_class = trtclass),
  set_agd_arm(pso_agd,
              study = studyc,
              trt = trtc,
              r = pasi75_r,
              n = pasi75_n,
              trt_class = trtclass)
)

# Print network details
pso_net

# Add integration points to the network
pso_net <- add_integration(pso_net,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
  weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
  psa = distr(qbern, prob = psa),
  n_int = 64)

\donttest{
# Fitting a ML-NMR model.
# Specify a regression model to include effect modifier interactions for five
# covariates, along with main (prognostic) effects. We use a probit link and
# specify that the two-parameter Binomial approximation for the aggregate-level
# likelihood should be used. We set treatment-covariate interactions to be equal
# within each class. We narrow the possible range for random initial values with
# init_r = 0.1, since probit models in particular are often hard to initialise.
# Using the QR decomposition greatly improves sampling efficiency here, as is
# often the case for regression models.
pso_fit <- nma(pso_net, 
               trt_effects = "fixed",
               link = "probit",
               likelihood = "bernoulli2",
               regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
               class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10),
               init_r = 0.1,
               QR = TRUE)
pso_fit
}

## Newly-diagnosed multiple myeloma NMA
# Set up newly-diagnosed multiple myeloma network

head(ndmm_ipd)
head(ndmm_agd)

ndmm_net <- combine_network(
  set_ipd(ndmm_ipd,
          study, trt,
          Surv = Surv(eventtime / 12, status)),
  set_agd_surv(ndmm_agd,
               study, trt,
               Surv = Surv(eventtime / 12, status),
               covariates = ndmm_agd_covs))
\donttest{
# Fit Weibull (PH) model
ndmm_fit <- nma(ndmm_net, 
                likelihood = "weibull",
                prior_intercept = normal(scale = 100),
                prior_trt = normal(scale = 10),
                prior_aux = half_normal(scale = 10))

ndmm_fit
}
}
\references{
\insertAllCited{}
}
