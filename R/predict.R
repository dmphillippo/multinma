#' Predictions of absolute effects from NMA models
#'
#' Obtain predictions of absolute effects from NMA models fitted with [nma()].
#' For example, if a model is fitted to binary data with a logit link, predicted
#' outcome probabilities or log odds can be produced.
#'
#' @param object A `stan_nma` object created by [nma()]
#' @param baseline An optional [distr()] distribution for the baseline response
#'   (i.e. intercept) on the linear predictor scale, about which to produce
#'   absolute effects. For example, in a model with a logit link, this would be
#'   a distribution for the baseline log odds of an event. If `NULL`,
#'   predictions are produced using the baseline response for each study in the
#'   network with IPD or contrast-based AgD.
#' @param newdata Only required if a regression model is fitted and `baseline`
#'   is specified. A data frame of covariate details, for which to produce
#'   predictions. Column names must match variables in the regression model.
#'
#'   If `type = "aggregate"` this should either be a data frame with integration
#'   points as produced by [add_integration()] (one row per study), or a data
#'   frame with individual covariate values (one row per individual) which are
#'   summarised over.
#'
#'   If `type = "individual"` this should be a data frame of individual
#'   covariate values, one row per individual.
#'
#'   If `NULL`, prections are produced for all studies with IPD and/or
#'   contrast-based AgD in the network, depending on the value of `type`.
#' @param study Column of `newdata` which specifies study names or IDs, required
#'   if `newdata` is specified.
#' @param type Whether to produce predictions on the `"link"` scale (the
#'   default, e.g. log odds) or `"response"` scale (e.g. probabilities).
#' @param level The level at which predictions are produced, either
#'   `"aggregate"` (the default), or `"individual"`. If `baseline` is not
#'   specified, predictions are produced for all IPD studies in the network if
#'   `type` is `"individual"` or `"aggregate"`, and for all contrast-based AgD
#'   studies in the network if `type` is `"aggregate"`.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @examples
predict.stan_nma <- function(object,
                             baseline = NULL, newdata = NULL, study = NULL,
                             type = c("link", "response"),
                             level = c("aggregate", "individual"),
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                             summary = TRUE) {

}
