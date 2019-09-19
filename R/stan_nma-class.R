#' The stan_nma class
#'
#' The `stan_nma` class contains the results from running a model with the
#' function [nma()].
#'
#' @rdname stan_nma-class
#' @name stan_nma-class
#' @aliases stan_nma
#'
#' @details Objects of class `stan_nma` have the following components:
#'   \describe{
#'   \item{`network`}{The network data from which the model was run (class
#'   [nma_data] or [mlnmr_data])}
#'   \item{`stanfit`}{The `stanfit` object returned by calling
#'   [rstan::sampling()] for the model}
#'   }
#'
NULL

#' Print stan_nma objects
#'
#' @param x A [stan_nma] object
#' @param ... Further arguments passed to [print.stanfit()]
#'
#' @export
print.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  dots <- list(...)
  include <- "pars" %in% names(dots)
  dots <- rlang::dots_list(x = sf,
                           pars = c("log_lik", "resdev", "theta_bar_cum", "theta2_bar_cum"),
                           include = include, !!! dots,
                           .homonyms = "last")
  do.call(print, dots)
  invisible(x)
}

#' as.stanfit
#'
#' Attempt to turn an object into a [rstan::stanfit] object.
#'
#' @param x an object
#' @param ... additional arguments
#'
#' @return A [rstan::stanfit] object.
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
