% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspline.R
\name{dmspline}
\alias{dmspline}
\alias{pmspline}
\alias{qmspline}
\alias{hmspline}
\alias{Hmspline}
\alias{rmst_mspline}
\title{Distribution functions for M-spline baseline hazards}
\usage{
dmspline(x, basis, scoef, rate, log = FALSE)

pmspline(q, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE)

qmspline(p, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE)

hmspline(x, basis, scoef, rate, log = FALSE)

Hmspline(x, basis, scoef, rate, log = FALSE)

rmst_mspline(t, basis, scoef, rate, start = 0)
}
\arguments{
\item{x, q}{Vector of quantiles}

\item{basis}{M-spline basis produced by \code{\link[splines2:mSpline]{splines2::mSpline()}}}

\item{scoef}{Vector (or matrix) of spline coefficients with length (or number
of columns) equal to the dimension of \code{basis}}

\item{rate}{Vector of rate parameters}

\item{log, log.p}{Logical; if \code{TRUE}, probabilities \code{p} are given as
\eqn{\log(p)}}

\item{lower.tail}{Logical; if \code{TRUE} (the default), probabilities are
\eqn{P(X \le x)}, otherwise \eqn{P(X > x)}}

\item{p}{Vector of probabilities}

\item{t}{Vector of times to which the restricted mean survival time is
calculated}

\item{start}{Optional left-truncation time or times. The returned restricted
mean survival will be conditioned on survival up to this time}
}
\value{
\code{dmspline()} gives the density, \code{pmspline()} gives the distribution
function (CDF), \code{qmspline()} gives the quantile function (inverse-CDF),
\code{hmspline()} gives the hazard function, \code{Hmspline()} gives the cumulative
hazard function, and \code{rmst_mspline()} gives restricted mean survival times.
}
\description{
Density, distribution, quantile, hazard, cumulative hazard, and restricted
mean survival time functions for the M-spline baseline hazards model.
}
\details{
Survival models with a flexible M-spline on the baseline hazard are
described by \insertCite{Brilleman2020;textual}{multinma}.
Piecewise-exponential baseline hazards are a special case where the degree
of the M-spline polynomial is 0.

The d/p/h/H functions are calculated from their definitions. \code{qmspline()}
uses numerical inversion via \code{\link[flexsurv:qgeneric]{flexsurv::qgeneric()}}. \code{rmst_mspline()}uses
numerical integration via \code{\link[flexsurv:rmst_generic]{flexsurv::rmst_generic()}}, except for the
special case of the piecewise-exponential hazard (i.e. degree 0 M-splines)
which uses the explicit formula from
\insertCite{Royston2013;textual}{multinma}.

Beyond the boundary knots, the hazard is assumed to be constant. (This
differs from the approach in \code{\link[splines2:mSpline]{splines2::mSpline()}} that extrapolates the
polynomial basis functions, which is numerically unstable and highly
dependent on the data just before the boundary knots.) As with all
extrapolation, care should be taken when evaluating the splines at times
beyond the boundary knots (either directly through the d/p/h/H/rmst
functions, or indirectly by requesting quantiles with \code{qmspline()} that
correspond to times beyond the boundary knots). For this reason evaluating
the (unrestricted) mean survival time is not generally recommended as this
requires integrating over an infinite time horizon (i.e. \code{rmst_mspline()}
with \code{t = Inf}).
}
\references{
\insertAllCited{}
}
