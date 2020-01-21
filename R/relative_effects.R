#' Relative treatment effects
#'
#' Generate (population-average) relative treatment efects. If a ML-NMR or
#' meta-regression model was fitted, these are specific to each study
#' population.
#'
#' @param x A `stan_nma` object created by [nma()]
#' @param newdata Only used if a regression model is fitted. A data frame of
#'   study details, one row per study, giving the covariate values at which to
#'   produce relative effects. Column names must match variables in the
#'   regression model. If `NULL`, relative effects are produced for all studies
#'   in the network.
#' @param study Column of `newdata` which specifies study names, otherwise
#'   studies will be labelled by row number.
#' @param all_contrasts Logical, generate estimates for all contrasts (`TRUE`),
#'   or just the "basic" contrasts against the network reference treatment
#'   (`FALSE`)? Default `FALSE`.
#'
#' @return A [nma_summary] object.
#' @export
#'
#' @examples
relative_effects <- function(x, newdata = NULL, study = NULL, all_contrasts = FALSE) {

}
