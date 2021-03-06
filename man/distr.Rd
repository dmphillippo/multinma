% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{distr}
\alias{distr}
\title{Specify a general marginal distribution}
\usage{
distr(qfun, ...)
}
\arguments{
\item{qfun}{an inverse CDF, either as a function name or a string}

\item{...}{parameters of the distribution as arguments to \code{qfun}, these will
be quoted and evaluated later in the context of the aggregate data sources}
}
\value{
An object of class \link{distr}.
}
\description{
\code{distr()} is used within the function \code{\link[=add_integration]{add_integration()}} to specify marginal
distributions for the covariates, via a corresponding inverse CDF. It is also
used in \code{\link[=predict.stan_nma]{predict.stan_nma()}} to specify a distribution for the baseline
response (intercept) when predicting absolute outcomes.
}
\details{
The function \code{qfun} should have a formal argument called \code{p}. This
restriction serves as a crude check for inverse CDFs (e.g. an error will be
given if \code{dnorm} is used instead of \code{qnorm}). If a user-written CDF is
supplied, it must have an argument \code{p} which takes a vector of
probabilities.
}
\examples{
## Specifying marginal distributions for integration

df <- data.frame(x1_mean = 2, x1_sd = 0.5, x2 = 0.8)

# Distribution parameters are evaluated in the context of the data frame
add_integration(df,
                x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                x2 = distr(qbern, prob = x2),
                cor = diag(2))

}
\seealso{
\code{\link[=add_integration]{add_integration()}} where \code{distr()} is used to specify marginal
distributions for covariates to integrate over, and \code{\link[=predict.stan_nma]{predict.stan_nma()}}
where \code{distr()} is used to specify a distribution on the baseline response.
}
