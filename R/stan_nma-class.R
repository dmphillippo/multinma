#' The stan_nma class
#'
#' The `stan_nma` and `stan_mlnmr` classes contains the results from running a
#' model with the function [nma()].
#'
#' @rdname stan_nma-class
#' @name stan_nma-class
#' @aliases stan_nma stan_mlnmr
#'
#' @details Objects of class `stan_nma` and `stan_mlnmr` have the following
#'   components:
#'   \describe{
#'   \item{`network`}{The network data from which the model was run (class
#'   [nma_data] for `stan_nma`, or class [mlnmr_data] for `stan_mlnmr`)}
#'   \item{`stanfit`}{The `stanfit` object returned by calling
#'   `\link[rstan:stanmodel-method-sampling]{sampling()}` for the model}
#'   \item{`trt_effects`}{Whether fixed or random effects were used (character
#'   string)}
#'   \item{`consistency`}{The consistency/inconsistency model used (character
#'   string)}
#'   \item{`regression`}{The regression model used (formula)}
#'   \item{`xbar`}{A named vector of values used for centering}
#'   \item{`likelihood`}{The likelihood used (character string)}
#'   \item{`link`}{The link function used (character string)}
#'   \item{`priors`}{A list containing the priors used (as [nma_prior] objects)}
#'   }
#'
#' The `stan_mlnmr` sub-class inherits from `stan_nma`, and differs only in the
#' class of the `network` object.
#'
NULL

#' Print stan_nma objects
#'
#' @param x A [stan_nma] object
#' @param ... Further arguments passed to [print.stanfit()]
#'
#' @export
print.stan_nma <- function(x, ...) {
  if (inherits(x$network, "mlnmr_data")) type <- "ML-NMR"
  else type <- "NMA"
  cglue("A {x$trt_effects} effects {type} with a {x$likelihood} likelihood ({x$link} link).")
  if (x$consistency != "consistency") cglue("An inconsistency model ('{x$consistency}') was fitted.")
  if (!is.null(x$regression)) cglue("Regression model: {rlang::as_label(x$regression)}.")

  sf <- as.stanfit(x)
  dots <- list(...)
  include <- "pars" %in% names(dots)
  dots <- rlang::dots_list(x = sf,
                           pars = c("log_lik", "resdev", "fitted",
                                    "theta_bar_cum", "theta2_bar_cum",
                                    "mu", "delta"),
                           include = include, !!! dots,
                           .homonyms = "last")
  do.call(print, dots)
  invisible(x)
}

#' as.stanfit
#'
#' Attempt to turn an object into a `\link[rstan:stanfit-class]{stanfit}` object.
#'
#' @param x an object
#' @param ... additional arguments
#'
#' @return A `\link[rstan:stanfit-class]{stanfit}` object.
#' @export
#'
#' @examples
as.stanfit <- function(x, ...) {
  UseMethod("as.stanfit")
}

#' @export
#' @noRd
as.stanfit.stan_nma <- function(x, ...) {
  return(x[["stanfit"]])
}

#' @export
#' @noRd
as.stanfit.default <- function(x, ...) {
  abort(glue::glue("Cannot coerce object of class '{class(x)}' to 'stanfit'."))
}

#' as.array
#'
#' Turn a `stan_nma` object into a 3D array \[Iteration, Chain, Parameter\].
#' Enables `\link[bayesplot:bayesplot-package]{bayesplot}` functions to
#' seamlessly work on `stan_nma` objects.
#'
#' @param x an object
#' @param ... additional arguments to [as.array.stanfit]
#'
#' @export
as.array.stan_nma <- function(x, ...) {
  return(as.array(as.stanfit(x), ...))
}

#' as.data.frame
#'
#' Turn a `stan_nma` object into a data frame containing posterior samples (post
#' warm-up) of each parameter.
#'
#' @param x an object
#' @param ... additional arguments to [as.data.frame.stanfit]
#'
#' @export
as.data.frame.stan_nma <- function(x, ...) {
  return(as.data.frame(as.stanfit(x), ...))
}

#' as.matrix
#'
#' Turn a `stan_nma` object into a matrix containing posterior samples (post
#' warm-up) of each parameter.
#'
#' @param x an object
#' @param ... additional arguments to [as.matrix.stanfit]
#'
#' @export
as.matrix.stan_nma <- function(x, ...) {
  return(as.matrix(as.stanfit(x), ...))
}

#' Model comparison using the `loo` package
#'
#' The `\link[loo:loo]{loo()}` and `\link[loo:waic]{waic()}` functions from the `loo`
#' package may be called directly on [stan_nma] and [stan_mlnmr] objects.
#'
#' @param x An object of class [stan_nma] or [stan_mlnmr]
#' @param ... Further arguments to `\link[rstan:stanfit-method-loo]{loo()}` or
#'   `\link[loo:waic]{waic()}`
#'
#' @export
#' @importFrom loo loo
#' @rdname loo
#' @aliases loo
loo.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  return(rstan::loo(sf, ...))
}

#' @export
#' @importFrom loo waic
#' @rdname loo
#' @aliases waic
waic.stan_nma <- function(x, ...) {
  ll <- as.array(x, pars = "log_lik")
  return(loo::waic(ll, ...))
}

#' Matrix of plots for a `stan_nma` object
#'
#' A [pairs()] method for `stan_nma` objects, which calls
#' `\link[rstan:stanfit-method-pairs]{pairs()}` on the underlying `stanfit`
#' object.
#'
#' @param x An object of class `stan_nma`
#' @param ... Other arguments passed to
#'   `\link[rstan:stanfit-method-pairs]{pairs()}`, such as `pars` to select the
#'   parameters to display.
#'
#' @return
#' @export
#'
#' @examples
pairs.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  pairs(sf, ...)
}
