#' multinma: A Package for Network Meta-Analysis of Individual and Aggregate
#' Data in Stan
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
#' @importFrom dplyr %>%
#' @importFrom rlang abort warn inform enquo .data :=
#' @importFrom rstan sampling
#' @importFrom Rdpack reprompt
#' @importFrom graphics plot pairs
#' @importFrom grDevices nclass.Sturges
#' @importFrom stats complete.cases sd median quantile model.frame model.matrix
#'   model.offset terms optim pbinom dbinom qbinom as.formula update.formula
#'   weighted.mean runif plogis pnorm
#' @importFrom utils packageVersion
#'
#' @references
#'
NULL

# Stop R CMD check thinking . used in pipes is an undeclared global variable
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))
