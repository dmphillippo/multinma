## Template argument docs for data functions set_*
#' @param data a data frame
#' @param trt column of \code{data} specifying treatments, coded using integers,
#'   strings, or factors
#' @param study column of \code{data} specifying the studies, coded using
#'   integers, strings, or factors
#' @param y column of \code{data} specifying a continuous outcome
#' @param se column of \code{data} specifying the standard error for a
#'   continuous outcome
#' @param r column of \code{data} specifying a binary or Binomial outcome count
#' @param n column of \code{data} specifying Binomial outcome numerator
#' @param E column of \code{data} specifying the total time at risk for Poisson
#'   outcomes
#' @param surv column of \code{data} specifying a survival outcome, using the
#'   \code{survival::Surv} function
