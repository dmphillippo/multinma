#' Network meta-analysis models
#'
#' The `nma` function fits network meta-analysis and (multilevel) network
#' meta-regression models in Stan.
#'
#' @param network An `nma_data` object, as created by the functions `set_*()`,
#'   `combine_network()`, or `add_integration()`
#' @param consistency Character string specifying the type of (in)consistency
#'   model to fit, either `"consistency"`, `"nodesplit"`, or `"ume"`
#' @param trt_effects Character string specifying either `"fixed"` or `"random"` effects
#' @param regression A one-sided model formula, specifying the prognostic and
#'   effect-modifying terms for a regression model
#' @param likelihood A suitable likelihood specification, if unspecified will be
#'   inferred from the data
#' @param ... Further arguments passed to [rstan::sampling()], such as `iter`,
#'   `chains`, `cores`, etc.
#' @param prior_intercept Specification of prior distribution for the intercept
#' @param prior_trt Specification of prior distribution for the treatment effects
#' @param prior_het Specification of prior distribution for the heterogeneity
#'   standard deviation (if `trt_effects = "random"`)
#' @param prior_reg Specification of prior distribution for the regression
#'   coefficients (if `regression` formula specified)
#' @param prior_aux Specification of prior distribution for the auxilliary
#'   parameter, if applicable
#' @param QR Logical scalar (default `FALSE`), whether to apply a QR
#'   decomposition to the model design matrix
#' @param adapt_delta See [adapt_delta] for details
#'
#' @return A [stan_nma] object.
#' @export
#'
#' @examples
nma <- function(network,
                consistency = c("consistency", "nodesplit", "ume"),
                trt_effects = c("fixed", "random"),
                regression = NULL,
                likelihood = NULL,
                ...,
                prior_intercept = normal(sd = 10),
                prior_trt = normal(sd = 10),
                prior_het = normal(sd = 5),
                prior_reg = normal(sd = 10),
                prior_aux = normal(sd = 5),
                QR = FALSE,
                adapt_delta = NULL) {

}


#' @param ipd_x Design matrix for IPD studies
#' @param ipd_y Outcome vector for IPD studies
#' @param agd_arm_x  Design matrix for AgD studies (arm-based)
#' @param agd_arm_y  Outcome vector for AgD studies (arm-based)
#' @param agd_contrast_x  Design matrix for AgD studies (contrast-based)
#' @param agd_contrast_y  Outcome vector for AgD studies (contrast-based)
#'
#' @export
#'
#' @rdname nma
nma.fit <- function(ipd_x, ipd_y,
                    agd_arm_x, agd_arm_y,
                    agd_contrast_x, agd_contrast_y,
                    consistency = c("consistency", "nodesplit", "ume"),
                    trt_effects = c("fixed", "random"),
                    regression = NULL,
                    likelihood = NULL,
                    ...,
                    prior_intercept = normal(sd = 10),
                    prior_trt = normal(sd = 10),
                    prior_het = normal(sd = 5),
                    prior_reg = normal(sd = 10),
                    prior_aux = normal(sd = 5),
                    QR = FALSE,
                    adapt_delta = NULL) {

}
