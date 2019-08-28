#' Add numerical integration points to aggregate data
#'
#' This function creates numerical integration points using a Gaussian copula
#' approach, as described in \insertCite{methods_paper}{multinma}.
#'
#' @param network an `nma_data` object, as created by the `set_*()` functions or
#'   `combine_network()`
#' @param ... distributions for covariates, see "Details"
#' @param cor correlation matrix to use for generating the integration points.
#'   By default, this takes a weighted correlation matrix from all IPD studies.
#' @param n_int number of integration points to generate, default 100
#' @param int_args a named list of arguments to pass to [randtoolbox::sobol()]
#'
#' @return An object of class [nma_data].
#' @export
#'
#' @details The arguments passed to `...` specify distributions for the
#'   covariates. Argument names specify the name of the covariate, which should
#'   match a covariate name in the IPD (if IPD are present). The required
#'   marginal distribution is then specified using the function [distr()].
#'
#' @examples
add_integration <- function(network, ..., cor = NULL, n_int = 100L, int_args = list()) {
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the `set_*` or `combine_network` functions.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  if (!is.numeric(n_int) || n_int != trunc(n_int) || n_int <= 0 || length(n_int) > 1) {
    abort("`n_int` should be a positive integer.")
  }

  if (!is.list(int_args) || (!rlang::is_empty(int_args) && !rlang::is_named(int_args))) {
    abort("`int_args` should be a named list.")
  }
}


#' Specify a general marginal distribution
#'
#' `distr()` is used within the function [add_integration()] to specify marginal
#' distributions for the covariates, via a corresponding inverse CDF.
#'
#' @param qfun an inverse CDF, either as a function name or a string
#' @param ... parameters of the distribution as arguments to `qfun`, these will
#'   be quoted and evaluated later in the context of the aggregate data sources
#'
#' @return An object of class [distr].
#' @export
#'
#' @examples
distr <- function(qfun, ...) {
  qfun <- match.fun(qfun)
  d <- as.list(rlang::enquos(...))
  d$qfun <- qfun
  class(d) <- "distr"
  return(d)
}
