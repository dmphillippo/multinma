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

  likelihood <- check_likelihood(likelihood, network$outcome)
  link <- check_link(link, likelihood)

  # Check priors
  if (!inherits(prior_intercept, "nma_prior")) abort("`prior_intercept` should be a prior distribution, see ?priors.")
  if (!inherits(prior_trt, "nma_prior")) abort("`prior_trt` should be a prior distribution, see ?priors.")
  if (!inherits(prior_het, "nma_prior")) abort("`prior_het` should be a prior distribution, see ?priors.")
  if (!inherits(prior_reg, "nma_prior")) abort("`prior_reg` should be a prior distribution, see ?priors.")
  if (!inherits(prior_aux, "nma_prior")) abort("`prior_aux` should be a prior distribution, see ?priors.")

  # Check other args
  if (!is.logical(QR) || length(QR) > 1) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!is.numeric(int_thin) ||
      length(int_thin) > 1 ||
      int_thin < 1 ||
      trunc(int_thin) != int_thin) abort("`int_thin` should be an integer > 0.")

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!is.numeric(adapt_delta) ||
      length(adapt_delta) > 1 ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Get design matrices and outcomes
  if (has_ipd(network)) {
    dat_ipd <- network$ipd
    o_ipd <- get_outcome_variables(dat_ipd, network$outcome$ipd)
  } else {
    dat_ipd <- tibble::tibble()
    o_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    dat_agd_arm <- network$agd_arm
    o_agd_arm <- get_outcome_variables(dat_agd_arm, network$outcome$agd_arm)

    # Set up integration variables if present
    if (inherits(network, "mlnmr_data")) {
      dat_agd_arm <- dat_agd_arm %>%
        dplyr::select(-dplyr::one_of(intersect(network$int_names, colnames(dat_agd_arm)))) %>%
        dplyr::rename_at(dplyr::vars(dplyr::starts_with(".int_")), ~gsub(".int_", "", ., fixed = TRUE)) %>%
        tidyr::unnest(!!! rlang::syms(network$int_names))
    }
  } else {
    dat_agd_arm <- tibble::tibble()
    o_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    dat_agd_contrast <- network$agd_contrast
    o_agd_contrast <- get_outcome_variables(dat_agd_contrast, network$outcome$agd_contrast)

    # Set up integration variables if present
    if (inherits(network, "mlnmr_data")) {
      dat_agd_contrast <- dat_agd_contrast %>%
        dplyr::select(-dplyr::one_of(union(network$int_names, colnames(dat_agd_contrast)))) %>%
        dplyr::rename_at(dplyr::vars(dplyr::starts_with(".int_")), ~gsub(".int_", "", ., fixed = TRUE)) %>%
        tidyr::unnest(!!! rlang::syms(network$int_names))
    }
  } else {
    dat_agd_contrast <- tibble::tibble()
    o_agd_contrast <- NULL
  }

  # Construct design matrix all together then split out, so that same dummy
  # coding is used everywhere
  dat_all <- dplyr::bind_rows(dat_ipd, dat_agd_arm, dat_agd_contrast)

  if (!is.null(regression)) {
    # Warn if combining AgD and IPD in a meta-regression without using integration
    if (!inherits(network, "mlnmr_data") && has_ipd(network) &&
        (has_agd_arm(network) || has_agd_contrast(network))) {
      warn(glue::glue("No integration points available, using naive plug-in model at aggregate level.\n",
                      "Use `add_integration()` to add integration points to the network."))
    }

    nma_formula <- update.formula(regression, ~-1 + .study + .trt + .)
  } else {
    nma_formula <- ~-1 + .study + .trt
  }

  X_all <- model.matrix(nma_formula, data = dat_all)

  if (has_ipd(network)) {
    X_ipd <- X_all[1:nrow(dat_ipd), ]
  } else {
    X_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    X_agd_arm <- X_all[nrow(dat_ipd) + 1:nrow(dat_agd_arm), ]
  } else {
    X_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    X_agd_contrast <- X_all[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast), ]

    # Need to difference .trtb terms - write custom model.matrix?
    abort("Contrast-based AgD not yet supported.")
  } else {
    X_agd_contrast <- NULL
  }


}


#' @param ipd_x Design matrix for IPD studies
#' @param ipd_y Outcome data frame for IPD studies
#' @param agd_arm_x  Design matrix for AgD studies (arm-based)
#' @param agd_arm_y  Outcome data frame for AgD studies (arm-based)
#' @param agd_contrast_x  Design matrix for AgD studies (contrast-based)
#' @param agd_contrast_y  Outcome data frame for AgD studies (contrast-based)
#' @param n_int Number of numerical integration points used
#' @param RE_cor Random effects correlation matrix, when `trt_effects = "random"`
#'
#' @export
#'
#' @rdname nma
nma.fit <- function(ipd_x, ipd_y,
                    agd_arm_x, agd_arm_y,
                    agd_contrast_x, agd_contrast_y,
                    n_int = NULL,
                    consistency = c("consistency", "nodesplit", "ume"),
                    trt_effects = c("fixed", "random"),
                    RE_cor = NULL,
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
  if (!is.numeric(int_thin) ||
      length(int_thin) > 1 ||
      int_thin < 1 ||
      trunc(int_thin) != int_thin) abort("`int_thin` should be an integer > 0.")

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!is.numeric(adapt_delta) ||
             length(adapt_delta) > 1 ||
             adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

}


check_likelihood <- function(x, outcome = NULL) {
  valid_lhood <- c("normal", "bernoulli", "binomial", "poisson")
  default_lhood <- list(binary = "bernoulli",
                        count = "binomial",
                        rate = "poisson",
                        continuous = "normal")

  # Choose default likelihood if NULL
  if (is.null(x) && is.null(outcome)) abort("Please specify a suitable likelihood.")
  else if (is.null(x)) {
    if (!is.na(outcome$ipd)) x <- default_lhood[[outcome$ipd]]
    else if (!is.na(outcome$agd_arm)) x <- default_lhood[[outcome$agd_arm]]
    else if (!is.na(outcome$agd_contrast)) x <- default_lhood[[outcome$agd_contrast]]
    else abort("No outcomes specified in network data.")

  # Check valid option if given
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_lhood) {
    abort(glue::glue("`likelihood` should be a character string specifying a valid likelihood.\n",
                     "Suitable options are currently: ",
                     glue::glue_collapse(dQuote(valid_lhood, FALSE), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

check_link <- function(x, lik) {
  valid_link <- list(normal = c("identity", "log"),
                     bernoulli = c("logit", "probit"),
                     binomial = c("logit", "probit"),
                     poisson = "log")[[lik]]

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

#' Get outcome variables from internal nma_data
#'
#' @param x A data frame
#' @param o_type Outcome type
#'
#' @return A data frame with outcome variables selected
#' @noRd
get_outcome_variables <- function(x, o_type) {
  o_vars <- list(
    binary = ".r",
    rate = c(".r", ".E"),
    count = c(".r", ".n"),
    continuous = c(".y", ".se")
  )[[o_type]]

  return(
    dplyr::select(x, dplyr::one_of(intersect(o_vars, colnames(x))))
  )
}
