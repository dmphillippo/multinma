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
