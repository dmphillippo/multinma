% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{add_integration}
\alias{add_integration}
\alias{add_integration.default}
\alias{add_integration.data.frame}
\alias{add_integration.nma_data}
\alias{unnest_integration}
\title{Add numerical integration points to aggregate data}
\usage{
add_integration(x, ...)

\method{add_integration}{default}(x, ...)

\method{add_integration}{data.frame}(x, ..., cor = NULL, n_int = 1000L, int_args = list())

\method{add_integration}{nma_data}(x, ..., cor = NULL, n_int = 1000L, int_args = list())

unnest_integration(data)
}
\arguments{
\item{x}{An \code{nma_data} object, as created by the \verb{set_*()} functions or
\code{combine_network()}, or data frame}

\item{...}{Distributions for covariates, see "Details"}

\item{cor}{Correlation matrix to use for generating the integration points.
By default, this takes a weighted correlation matrix from all IPD studies.
Rows and columns should match the order of covariates specified in \code{...}.}

\item{n_int}{Number of integration points to generate, default 1000}

\item{int_args}{A named list of arguments to pass to
\code{\link[randtoolbox:quasiRNG]{sobol()}}}

\item{data}{Data frame with nested integration points, stored in list
columns as \verb{.int_<variable name>}}
}
\value{
For the \code{nma_data} method, an object of class \link{nma_data}. For the
\code{data.frame} method, the input data frame is returned (as a \link{tibble}) with
an added column for each covariate (prefixed with ".int_"), containing the
numerical integration points nested as length-\code{n_int} vectors within each
row. For \code{unnest_integration()}, a data frame with integration points
unnested.
}
\description{
The \code{add_integration()} generic creates numerical integration points using a
Gaussian copula approach, as described in
\insertCite{methods_paper;textual}{multinma}. Methods are available for
networks stored in \code{nma_data} objects, and for data frames. The function
\code{unnest_integration()} unnests integration points stored in a data frame, to
aid plotting or other exploration.
}
\details{
The arguments passed to \code{...} specify distributions for the
covariates. Argument names specify the name of the covariate, which should
match a covariate name in the IPD (if IPD are present). The required
marginal distribution is then specified using the function \code{\link[=distr]{distr()}}.
}
\examples{
## Plaque psoriasis ML-NMR - network setup and adding integration points
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
  n_int = 1000)


## Adding integration points to a data frame, e.g. for prediction
# Define a data frame of covariate summaries
new_agd_int <- data.frame(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1)

# Adding integration points, using the weighted average correlation matrix
# computed for the plaque psoriasis network
new_agd_int <- add_integration(new_agd_int,
  durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
  prevsys = distr(qbern, prob = prevsys),
  bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
  weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
  psa = distr(qbern, prob = psa),
  cor = pso_net$int_cor,
  n_int = 1000)

new_agd_int

}
\references{
\insertAllCited{}
}
