#' The 'multinma' package.
#'
#' @description This package implements network meta-analysis (NMA) and network
#'   meta-regression (NMR) models for aggregate data (AgD), individual patient
#'   data (IPD), and mixtures of both IPD and AgD using multilevel NMR (ML-NMR).
#'   Models are estimated in a Bayesian framwork using Stan.
#'
#' @docType package
#' @name multinma-package
#' @aliases multinma
#' @useDynLib multinma, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references Stan Development Team (2019). RStan: the R interface to Stan. R
#' package version 2.19.2. https://mc-stan.org
#'
NULL
