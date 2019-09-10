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
#' @param likelihood Character string specifying a likelihood, if unspecified
#'   will be inferred from the data
#' @param link Character string specifying a link function, if unspecified will
#'   default to the canonical link
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
#' @param int_thin A single integer value, the thinning factor for returning
#'   cumulative estimates of integration error
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
                link = NULL,
                ...,
                prior_intercept = normal(sd = 10),
                prior_trt = normal(sd = 10),
                prior_het = half_normal(sd = 5),
                prior_reg = normal(sd = 10),
                prior_aux = normal(sd = 5),
                QR = FALSE,
                adapt_delta = NULL,
                int_thin = 100L) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # Check model arguments
  consistency <- rlang::arg_match(consistency)
  if (length(consistency) > 1) abort("`consistency` must be a single string.")
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")

  if (!is.null(regression) && !rlang::is_formula(regression, lhs = FALSE)) {
    abort("`regression` should be a one-sided formula.")
  }

  likelihood <- check_likelihood(likelihood)
  link <- check_link(link, likelihood)

  # Check priors
  if (!inherits(prior_intercept, "nma_prior")) abort("`prior_intercept` should be a prior distribution, see ?priors.")
  if (!inherits(prior_trt, "nma_prior")) abort("`prior_trt` should be a prior distribution, see ?priors.")
  if (!inherits(prior_het, "nma_prior")) abort("`prior_het` should be a prior distribution, see ?priors.")
  if (!inherits(prior_reg, "nma_prior")) abort("`prior_reg` should be a prior distribution, see ?priors.")
  if (!inherits(prior_aux, "nma_prior")) abort("`prior_aux` should be a prior distribution, see ?priors.")

  # Check other args
  if (!is.logical(QR) || length(QR) > 1) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!is.numeric(adapt_delta) ||
      length(adapt_delta) > 1 ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

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
                    link = NULL,
                    ...,
                    prior_intercept = normal(sd = 10),
                    prior_trt = normal(sd = 10),
                    prior_het = half_normal(sd = 5),
                    prior_reg = normal(sd = 10),
                    prior_aux = normal(sd = 5),
                    QR = FALSE,
                    adapt_delta = NULL,
                    int_thin = 100L) {

  # Check model arguments
  consistency <- rlang::arg_match(consistency)
  if (length(consistency) > 1) abort("`consistency` must be a single string.")
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")

  if (!is.null(regression) && !rlang::is_formula(regression, lhs = FALSE)) {
    abort("`regression` should be a one-sided formula.")
  }

  likelihood <- check_likelihood(likelihood)
  link <- check_link(link, likelihood)

  # Check priors
  if (!inherits(prior_intercept, "nma_prior")) abort("`prior_intercept` should be a prior distribution, see ?priors.")
  if (!inherits(prior_trt, "nma_prior")) abort("`prior_trt` should be a prior distribution, see ?priors.")
  if (!inherits(prior_het, "nma_prior")) abort("`prior_het` should be a prior distribution, see ?priors.")
  if (!inherits(prior_reg, "nma_prior")) abort("`prior_reg` should be a prior distribution, see ?priors.")
  if (!inherits(prior_aux, "nma_prior")) abort("`prior_aux` should be a prior distribution, see ?priors.")

  # Check other args
  if (!is.logical(QR) || length(QR) > 1) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!is.numeric(adapt_delta) ||
      length(adapt_delta) > 1 ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

}


check_likelihood <- function(x) {
  valid_lhood <- c("normal", "bernoulli", "binomial")

  if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_lhood) {
    abort(glue::glue("`likelihood` should be a character string specifying a valid likelihood.\n",
                     "Suitable options are currently: ",
                     glue::glue_collapse(dQuote(valid_lhood, FALSE), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

check_link <- function(x, lik) {
  valid_link <- list(normal = c("identity", "log"),
                     bernoulli = c("logit", "probit", "cloglog"),
                     binomial = c("logit", "probit", "cloglog"))[[lik]]

  if (is.null(x)) {
    x <- valid_link[1]
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_link) {
    abort(glue::glue("`link` should be a character string specifying a valid link function.\n",
                     "Suitable options for a {lik} likelihood are currently: ",
                     glue::glue_collapse(dQuote(valid_link, FALSE), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}
