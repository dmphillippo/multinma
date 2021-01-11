#' Network meta-analysis models
#'
#' The `nma` function fits network meta-analysis and (multilevel) network
#' meta-regression models in Stan.
#'
#' @param network An `nma_data` object, as created by the functions `set_*()`,
#'   `combine_network()`, or `add_integration()`
#' @param consistency Character string specifying the type of (in)consistency
#'   model to fit, currently either `"consistency"` or `"ume"`
#' @param trt_effects Character string specifying either `"fixed"` or `"random"` effects
#' @param regression A one-sided model formula, specifying the prognostic and
#'   effect-modifying terms for a regression model. Any references to treatment
#'   should use the `.trt` special variable, for example specifying effect
#'   modifier interactions as `variable:.trt` (see details).
#' @param class_interactions Character string specifying whether effect modifier
#'   interactions are specified as `"common"`, `"exchangeable"`, or
#'   `"independent"`.
#' @param likelihood Character string specifying a likelihood, if unspecified
#'   will be inferred from the data (see details)
#' @param link Character string specifying a link function, if unspecified will
#'   default to the canonical link (see details)
#' @param ... Further arguments passed to
#'   \code{\link[rstan:stanmodel-method-sampling]{sampling()}}, such as `iter`,
#'   `chains`, `cores`, etc.
#' @param prior_intercept Specification of prior distribution for the intercept
#' @param prior_trt Specification of prior distribution for the treatment effects
#' @param prior_het Specification of prior distribution for the heterogeneity
#'   (if `trt_effects = "random"`)
#' @param prior_het_type Character string specifying whether the prior
#'   distribution `prior_het` is placed on the heterogeneity standard deviation
#'   \eqn{\tau} (`"sd"`, the default), variance \eqn{\tau^2} (`"var"`), or
#'   precision \eqn{1/\tau^2} (`"prec"`).
#' @param prior_reg Specification of prior distribution for the regression
#'   coefficients (if `regression` formula specified)
#' @param prior_aux Specification of prior distribution for the auxiliary
#'   parameter, if applicable (see details)
#' @param QR Logical scalar (default `FALSE`), whether to apply a QR
#'   decomposition to the model design matrix
#' @param center Logical scalar (default `TRUE`), whether to center the
#'   (numeric) regression terms about the overall means
#' @param adapt_delta See [adapt_delta] for details
#' @param int_thin A single integer value, the thinning factor for returning
#'   cumulative estimates of integration error
#'
#' @details When specifying a model formula in the `regression` argument, the
#'   usual formula syntax is available (as interpreted by [model.matrix()]). The
#'   only additional requirement here is that the special variable `.trt` should
#'   be used to refer to treatment. For example, effect modifier interactions
#'   should be specified as `variable:.trt`. Prognostic (main) effects and
#'   interactions can be included together compactly as `variable*.trt`, which
#'   expands to `variable + variable:.trt` (plus `.trt`, which is already in the
#'   NMA model).
#'
#'   For the advanced user, the additional specials `.study` and `.trtclass` are
#'   also available, and refer to studies and (if specified) treatment classes
#'   respectively.
#'
#'   See \code{\link[multinma:priors]{?priors}} for details on prior
#'   specification. Default prior distributions are available, but may not be
#'   appropriate for the particular setting and will raise a warning if used. No
#'   attempt is made to tailor these defaults to the data provided. Please
#'   consider appropriate prior distributions for the particular setting,
#'   accounting for the scales of outcomes and covariates, etc. The function
#'   [plot_prior_posterior()] may be useful in examining the influence of the
#'   chosen prior distributions on the posterior distributions, and the
#'   \code{\link[multinma:summary.nma_prior]{summary()}} method for `nma_prior`
#'   objects prints prior intervals.
#'
#' @section Likelihoods and link functions:
#'   Currently, the following likelihoods and link functions are supported for
#'   each data type:
#'
#'   | \strong{Data type} | \strong{Likelihood}   | \strong{Link function} |
#'   |--------------------|-----------------------|------------------------|
#'   | \strong{Binary}    | `bernoulli`, `bernoulli2`| `logit`, `probit`, `cloglog`
#'   | \strong{Count}     | `binomial`, `binomial2`  | `logit`, `probit`, `cloglog`
#'   | \strong{Rate}      | `poisson`    | `log`
#'   | \strong{Continuous}| `normal`     | `identity`, `log`
#'   | \strong{Ordered}   | `ordered`    | `logit`, `probit`, `cloglog`
#'
#'   The `bernoulli2` and `binomial2` likelihoods correspond to a two-parameter
#'   Binomial likelihood for arm-based AgD, which more closely matches the
#'   underlying Poisson Binomial distribution for the summarised aggregate
#'   outcomes in a ML-NMR model than the typical (one parameter) Binomial
#'   distribution \insertCite{@see @methods_paper}{multinma}.
#'
#'   When a `cloglog` link is used, including an offset for log follow-up time
#'   (i.e. `regression = ~offset(log(time))`) results in a model on the log
#'   hazard \insertCite{@see @TSD2}{multinma}.
#'
#'   Further details on each likelihood and link function are given by
#'   \insertCite{TSD2;textual}{multinma}.
#'
#'
#' @section Auxiliary parameters:
#'   Auxiliary parameters are only present in the following models.
#'
#'   ## Normal likelihood with IPD
#'   When a Normal likelihood is fitted to IPD, the auxiliary parameters are the
#'   arm-level standard deviations \eqn{\sigma_{jk}} on treatment \eqn{k} in
#'   study \eqn{j}.
#'
#'   ## Ordered multinomial likelihood
#'   When fitting a model to \eqn{M} ordered outcomes, the auxiliary parameters
#'   are the latent cutoffs between each category, \eqn{c_0 < c_1 < \dots <
#'   c_M}. Only \eqn{c_2} to \eqn{c_{M-1}} are estimated; we fix \eqn{c_0 =
#'   -\infty}, \eqn{c_1 = 0}, and \eqn{c_M = \infty}. When specifying priors for
#'   these latent cutoffs, we choose to specify priors on the *differences*
#'   \eqn{c_{m+1} - c_m}. Stan automatically truncates any priors so that the
#'   ordering constraints are satisfied.
#'
#' @return `nma()` returns a [stan_nma] object, `nma.fit()` returns a [stanfit]
#'   object.
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' ## Smoking cessation NMA
#' @template ex_smoking_network
#' @template ex_smoking_nma_fe
#' @template ex_smoking_nma_re
#' @template ex_smoking_nma_re_ume
#' @examples
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_network
#' @template ex_plaque_psoriasis_integration
#' @template ex_plaque_psoriasis_mlnmr
#'
nma <- function(network,
                consistency = c("consistency", "ume"),
                trt_effects = c("fixed", "random"),
                regression = NULL,
                class_interactions = c("common", "exchangeable", "independent"),
                likelihood = NULL,
                link = NULL,
                ...,
                prior_intercept = .default(normal(scale = 100)),
                prior_trt = .default(normal(scale = 10)),
                prior_het = .default(half_normal(scale = 5)),
                prior_het_type = c("sd", "var", "prec"),
                prior_reg = .default(normal(scale = 10)),
                prior_aux = .default(),
                QR = FALSE,
                center = TRUE,
                adapt_delta = NULL,
                int_thin = max(network$n_int %/% 10, 1)) {

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

  if (is.null(network$classes)) {
    if (!missing(class_interactions)) {
      abort(paste("Setting `class_interactions` requires treatment classes to be specified in the network.",
                  "See set_*() argument `trt_class`.", sep = "\n"))
    } else if (!is.null(regression)) {
      inform(paste("Note: No treatment classes specified in network, any interactions in `regression` formula will be separate (independent) for each treatment.",
                   "Use set_*() argument `trt_class` and nma() argument `class_interactions` to change this.", sep = "\n"))
    }
  }
  class_interactions <- rlang::arg_match(class_interactions)
  if (length(class_interactions) > 1) abort("`class_interactions` must be a single string.")

  likelihood <- check_likelihood(likelihood, network$outcome)
  link <- check_link(link, likelihood)

  # When are priors on auxiliary parameters required?
  has_aux <- (likelihood == "normal" && has_ipd(network)) || likelihood == "ordered"

  # Are study intercepts present? Not if only contrast data
  has_intercepts <- has_agd_arm(network) || has_ipd(network)

  # Check priors
  if (!inherits(prior_intercept, "nma_prior")) abort("`prior_intercept` should be a prior distribution, see ?priors.")
  if (!inherits(prior_trt, "nma_prior")) abort("`prior_trt` should be a prior distribution, see ?priors.")
  if (!inherits(prior_het, "nma_prior")) abort("`prior_het` should be a prior distribution, see ?priors.")
  if (!inherits(prior_reg, "nma_prior")) abort("`prior_reg` should be a prior distribution, see ?priors.")
  if (!.is_default(prior_aux) && !inherits(prior_aux, "nma_prior")) abort("`prior_aux` should be a prior distribution, see ?priors.")

  prior_het_type <- rlang::arg_match(prior_het_type)


  # Prior defaults
  prior_defaults <- list()
  if (has_intercepts && .is_default(prior_intercept))
    prior_defaults$prior_intercept <- get_prior_call(prior_intercept)
  if (.is_default(prior_trt))
    prior_defaults$prior_trt <- get_prior_call(prior_trt)
  if (trt_effects == "random" && .is_default(prior_het))
    prior_defaults$prior_het <- get_prior_call(prior_het)
  if (!is.null(regression) && !is_only_offset(regression) && .is_default(prior_reg))
    prior_defaults$prior_reg <- get_prior_call(prior_reg)
  if (has_aux && .is_default(prior_aux)) {
    if (likelihood == "normal" && has_ipd(network)) {
      prior_aux <- .default(half_normal(scale = 5))
    } else if (likelihood == "ordered") {
      prior_aux <- .default(flat())
    }
    prior_defaults$prior_aux <- get_prior_call(prior_aux)
  }

  # Warn where default priors are used
  if (!rlang::is_empty(prior_defaults)) {
    warn(glue::glue(
      "Prior distributions were left at default values:",
      paste(paste(names(prior_defaults), prior_defaults, sep = " = "), collapse = "\n"),
      .sep = "\n"
    ))
  }

  # Check other args
  if (!rlang::is_bool(QR)) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_bool(center)) abort("`center` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_scalar_integerish(int_thin) ||
      int_thin < 1) abort("`int_thin` should be an integer >= 1.")

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!rlang::is_scalar_double(adapt_delta) ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Use numerical integration? TRUE if class mlnmr_data and regression is not NULL
  # (Avoids unnecessary use of integration points if regression formula not specified)
  use_int <- inherits(network, "mlnmr_data") && !is.null(regression)

  # Number of numerical integration points
  # Set to 1 if no numerical integration, so that regression on summary data is possible
  n_int <- if (use_int) network$n_int else 1

  # Warn if combining AgD and IPD in a meta-regression without using integration
  if (!is.null(regression) && !inherits(network, "mlnmr_data") && has_ipd(network) &&
      (has_agd_arm(network) || has_agd_contrast(network))) {
    warn(glue::glue("No integration points available, using naive plug-in model at aggregate level.\n",
                    "Use `add_integration()` to add integration points to the network."))
  }

  # Notify if default reference treatment is used
  if (.is_default(network$treatments))
    inform(glue::glue('Note: Setting "{levels(network$treatments)[1]}" as the network reference treatment.'))

  # Notify if network is disconnected
  if (!is_network_connected(network))
    inform("Note: Network is disconnected. See ?is_network_connected for more details.")

  # Get data for design matrices and outcomes
  if (has_ipd(network)) {
    dat_ipd <- network$ipd

    # Only take necessary columns
    dat_ipd <- get_model_data_columns(dat_ipd,
                                      regression = regression,
                                      label = "IPD")

    y_ipd <- get_outcome_variables(dat_ipd, network$outcome$ipd)
  } else {
    dat_ipd <- tibble::tibble()
    y_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    dat_agd_arm <- network$agd_arm

    y_agd_arm <- get_outcome_variables(dat_agd_arm, network$outcome$agd_arm)

    # Set up integration variables if present
    if (use_int) {
      idat_agd_arm <- .unnest_integration(dat_agd_arm)
    } else {
      idat_agd_arm <- dat_agd_arm
    }

    # Only take necessary columns
    idat_agd_arm <- get_model_data_columns(idat_agd_arm,
                                           regression = regression,
                                           label = "AgD (arm-based)")

  } else {
    dat_agd_arm <- idat_agd_arm <- tibble::tibble()
    y_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    dat_agd_contrast <- network$agd_contrast

    y_agd_contrast <- get_outcome_variables(dat_agd_contrast, network$outcome$agd_contrast)

    # Set up integration variables if present
    if (use_int) {
      idat_agd_contrast <-  .unnest_integration(dat_agd_contrast)
    } else {
      idat_agd_contrast <- dat_agd_contrast
    }

    # Only take necessary columns
    idat_agd_contrast <- get_model_data_columns(idat_agd_contrast,
                                                regression = regression,
                                                label = "AgD (contrast-based)")

    # Get covariance structure for relative effects, using .se on baseline arm
    Sigma_agd_contrast <- make_Sigma(dat_agd_contrast)

    # Split into baseline and non-baseline arms
    dat_agd_contrast_bl <- dplyr::filter(dat_agd_contrast, is.na(.data$.y))
    idat_agd_contrast_bl <- dplyr::filter(idat_agd_contrast, is.na(.data$.y))
    dat_agd_contrast_nonbl <- dplyr::filter(dat_agd_contrast, !is.na(.data$.y))
    idat_agd_contrast_nonbl <- dplyr::filter(idat_agd_contrast, !is.na(.data$.y))
    y_agd_contrast <- dplyr::filter(y_agd_contrast, !is.na(.data$.y))
  } else {
    dat_agd_contrast <- idat_agd_contrast <-
      dat_agd_contrast_bl <- idat_agd_contrast_bl <-
      dat_agd_contrast_nonbl <- idat_agd_contrast_nonbl <- tibble::tibble()
    y_agd_contrast <- NULL
    Sigma_agd_contrast <- NULL
  }

  # Combine
  idat_all <- dplyr::bind_rows(dat_ipd, idat_agd_arm, idat_agd_contrast_nonbl)
  idat_all_plus_bl <- dplyr::bind_rows(dat_ipd, idat_agd_arm, idat_agd_contrast)

  # Get sample sizes for centering
  if (!is.null(regression) && !is_only_offset(regression) && center) {
    # Check that required variables are present in each data set, and non-missing
    check_regression_data(regression,
                          dat_ipd = dat_ipd,
                          dat_agd_arm = idat_agd_arm,
                          dat_agd_contrast = idat_agd_contrast)

    # If IPD or IPD+AgD use weighted means for centering, otherwise with only AgD use raw mean

    if (has_ipd(network) && (has_agd_arm(network) || has_agd_contrast(network)) && !has_agd_sample_size(network))
      abort(paste("AgD study sample sizes not specified in network, cannot calculate centering values.",
                  "Specify `sample_size` in set_agd_*(), or set center = FALSE.", sep = "\n"))

    if (has_agd_arm(network)) {
      if (has_ipd(network)) {
        N_agd_arm <- network$agd_arm[[".sample_size"]]
      } else {
        N_agd_arm <- rep(n_int, nrow(network$agd_arm))
      }
    } else {
      N_agd_arm <- NULL
    }

    if (has_agd_contrast(network)) {
      if (has_ipd(network)) {
        N_agd_contrast <- network$agd_contrast[[".sample_size"]]
      } else {
        N_agd_contrast <- rep(n_int, nrow(network$agd_contrast))
      }
    } else {
      N_agd_contrast <- NULL
    }

    # Apply weights across integration points
    wts <- c(rep(1, nrow(dat_ipd)),
             rep(N_agd_arm / n_int, each = n_int),
             rep(N_agd_contrast / n_int, each = n_int))

    # Center numeric columns used in regression model
    reg_names <- all.vars(regression)

    # Ignore any variable(s) used as offset(s)
    reg_terms <- terms(regression)

    if (length(attr(reg_terms, "offset"))) {
      off_terms <- rownames(attr(reg_terms, "factors"))[attr(reg_terms, "offset")]
      off_names <- all.vars(as.formula(paste("~", off_terms, sep = "+")))
      reg_names <- setdiff(reg_names, off_names)
    }

    reg_numeric <- purrr::map_lgl(idat_all[, reg_names], is.numeric)

    # Take weighted mean of all rows (including baseline rows for contrast data)
    if (any(reg_numeric)) {
      xbar <- purrr::map_dbl(idat_all_plus_bl[, reg_names[reg_numeric]], weighted.mean, w = wts)
    } else {
      xbar <- NULL
    }

  } else {
    xbar <- NULL
  }

  # Make NMA formula
  nma_formula <- make_nma_formula(regression,
                                  consistency = consistency,
                                  classes = !is.null(network$classes),
                                  class_interactions = class_interactions)

  # Construct model matrix
  X_list <- make_nma_model_matrix(nma_formula = nma_formula,
                                  dat_ipd = dat_ipd,
                                  dat_agd_arm = idat_agd_arm,
                                  dat_agd_contrast = idat_agd_contrast,
                                  agd_contrast_bl = is.na(idat_agd_contrast$.y),
                                  xbar = xbar,
                                  consistency = consistency,
                                  classes = !is.null(network$classes))

  X_ipd <- X_list$X_ipd
  X_agd_arm <- X_list$X_agd_arm
  X_agd_contrast <- X_list$X_agd_contrast

  offset_ipd <- X_list$offset_ipd
  offset_agd_arm <- X_list$offset_agd_arm
  offset_agd_contrast <- X_list$offset_agd_contrast

  # Construct RE correlation matrix
  if (trt_effects == "random") {

    # Get study/treatment data
    if (has_ipd(network)) {
      tdat_ipd_arm <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt)
    } else {
      tdat_ipd_arm <- tibble::tibble()
    }

    if (has_agd_arm(network)) {
      tdat_agd_arm <- dplyr::select(dat_agd_arm, .data$.study, .data$.trt)
    } else {
      tdat_agd_arm <- tibble::tibble()
    }

    if (has_agd_contrast(network)) {
      tdat_agd_contrast_nonbl <- dplyr::select(dat_agd_contrast_nonbl, .data$.study, .data$.trt)
    } else {
      tdat_agd_contrast_nonbl <- tibble::tibble()
    }

    tdat_all <- dplyr::bind_rows(tdat_ipd_arm, tdat_agd_arm, tdat_agd_contrast_nonbl)

    contr <- rep(c(FALSE, FALSE, TRUE),
                 times = c(nrow(tdat_ipd_arm), nrow(tdat_agd_arm), nrow(tdat_agd_contrast_nonbl)))

    if (consistency == "consistency") {
      .RE_cor <- RE_cor(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "reftrt")
      .which_RE <- which_RE(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "reftrt")
    } else if (consistency == "ume") {
      .RE_cor <- RE_cor(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "blshift")
      .which_RE <- which_RE(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "blshift")
    } else {
      abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
    }
  } else {
    .RE_cor <- NULL
    .which_RE <- NULL
  }

  # Fit using nma.fit
  stanfit <- nma.fit(ipd_x = X_ipd, ipd_y = y_ipd,
    agd_arm_x = X_agd_arm, agd_arm_y = y_agd_arm,
    agd_contrast_x = X_agd_contrast, agd_contrast_y = y_agd_contrast,
    agd_contrast_Sigma = Sigma_agd_contrast,
    n_int = n_int,
    ipd_offset = offset_ipd,
    agd_arm_offset = offset_agd_arm,
    agd_contrast_offset = offset_agd_contrast,
    trt_effects = trt_effects,
    RE_cor = .RE_cor,
    which_RE = .which_RE,
    likelihood = likelihood,
    link = link,
    ...,
    prior_intercept = prior_intercept,
    prior_trt = prior_trt,
    prior_het = prior_het,
    prior_het_type = prior_het_type,
    prior_reg = prior_reg,
    prior_aux = prior_aux,
    QR = QR,
    adapt_delta = adapt_delta,
    int_thin = int_thin)

  # Make readable parameter names for generated quantities
  fnames_oi <- stanfit@sim$fnames_oi

  # Labels for fitted values
  ipd_data_labels <- if (has_ipd(network)) make_data_labels(dat_ipd$.study, dat_ipd$.trt) else NULL
  agd_arm_data_labels <- if (has_agd_arm(network)) make_data_labels(dat_agd_arm$.study, dat_agd_arm$.trt) else NULL

  if (has_agd_contrast(network)) {
    dat_agd_contrast_nonbl <-
      dplyr::left_join(dat_agd_contrast_nonbl,
                       dplyr::transmute(dat_agd_contrast_bl, .data$.study, .trt_b = .data$.trt),
                       by = ".study")

    agd_contrast_data_labels <- make_data_labels(dat_agd_contrast_nonbl$.study,
                                                 dat_agd_contrast_nonbl$.trt,
                                                 dat_agd_contrast_nonbl$.trt_b)
  } else {
    agd_contrast_data_labels <- NULL
  }

  if (has_ipd(network))
    fnames_oi[grepl("^fitted_ipd\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_ipd[", ipd_data_labels, "]")

  if (has_agd_arm(network))
    fnames_oi[grepl("^fitted_agd_arm\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_agd_arm[", agd_arm_data_labels, "]")

  if (has_agd_contrast(network))
    fnames_oi[grepl("^fitted_agd_contrast\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_agd_contrast[", agd_contrast_data_labels, "]")

  # Labels for RE deltas
  if (trt_effects == "random") {

    if (consistency == "ume") { # Baseline shift models (currently only UME)
      if (has_ipd(network)) {
        ipd_arms <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt) %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::mutate(.trt_b = sort(.data$.trt)[1])
        ipd_delta_labels <- make_data_labels(ipd_arms$.study, ipd_arms$.trt, ipd_arms$.trt_b)
      } else {
        ipd_delta_labels <- NULL
      }

      if (has_agd_arm(network)) {
        agd_arm_trt_b <- dat_agd_arm %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::mutate(.trt_b = sort(.data$.trt)[1]) %>%
          dplyr::pull(.data$.trt_b)
        agd_arm_delta_labels <- make_data_labels(dat_agd_arm$.study, dat_agd_arm$.trt, agd_arm_trt_b)
      } else {
        agd_arm_delta_labels <- NULL
      }

    } else { # Reference treatment models
      if (has_ipd(network)) {
        ipd_arms <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt)
        ipd_delta_labels <- make_data_labels(ipd_arms$.study, ipd_arms$.trt)
      } else {
        ipd_delta_labels <- NULL
      }
      agd_arm_delta_labels <- agd_arm_data_labels
    }

    agd_contrast_delta_labels <- agd_contrast_data_labels

    delta_labels <- c(ipd_delta_labels,
                      agd_arm_delta_labels,
                      agd_contrast_delta_labels)[.which_RE > 0]

    fnames_oi[grepl("^delta\\[[0-9]+\\]$", fnames_oi)] <- paste0("delta[", delta_labels, "]")
  }

  # Labels for log_lik, resdev (only one entry per AgD contrast study)
  dev_labels <- c(ipd_data_labels,
                  agd_arm_data_labels,
                  if (has_agd_contrast(network)) dat_agd_contrast_bl$.study else NULL)

  fnames_oi[grepl("^log_lik\\[[0-9]+\\]$", fnames_oi)] <- paste0("log_lik[", dev_labels, "]")
  fnames_oi[grepl("^resdev\\[[0-9]+\\]$", fnames_oi)] <- paste0("resdev[", dev_labels, "]")

  # Labels for cumulative integration points
  if (inherits(network, "mlnmr_data") && (has_agd_arm(network) || has_agd_contrast(network))) {
    n_int_thin <- n_int %/% int_thin

    if (has_agd_arm(network)) {
      agd_arm_cumint_labels <- rep(agd_arm_data_labels, each = n_int_thin)
      agd_arm_cumint_labels <- paste0(agd_arm_cumint_labels, ", ", rep_len(1:n_int_thin * int_thin, length.out = length(agd_arm_cumint_labels)))

      fnames_oi[grepl("^theta_bar_cum_agd_arm\\[[0-9]+\\]$", fnames_oi)] <- paste0("theta_bar_cum_agd_arm[", agd_arm_cumint_labels, "]")
      if (likelihood %in% c("bernoulli2", "binomial2"))
        fnames_oi[grepl("^theta2_bar_cum\\[[0-9]+\\]$", fnames_oi)] <- paste0("theta2_bar_cum[", agd_arm_cumint_labels, "]")
    }

    if (has_agd_contrast(network)) {
      agd_contrast_cumint_labels <- rep(agd_contrast_data_labels, each = n_int_thin)
      agd_contrast_cumint_labels <- paste0(agd_contrast_cumint_labels, ", ", rep_len(1:n_int_thin * int_thin, length.out = length(agd_contrast_cumint_labels)))

      fnames_oi[grepl("^theta_bar_cum_agd_contrast\\[[0-9]+\\]$", fnames_oi)] <- paste0("theta_bar_cum_agd_contrast[", agd_contrast_cumint_labels, "]")
    }

  }

  stanfit@sim$fnames_oi <- fnames_oi

  # Create stan_nma object
  out <- list(network = network,
              stanfit = stanfit,
              trt_effects = trt_effects,
              consistency = consistency,
              regression = regression,
              class_interactions = if (!is.null(regression) && !is.null(network$classes)) class_interactions else NULL,
              xbar = xbar,
              likelihood = likelihood,
              link = link,
              priors = list(prior_intercept = if (has_intercepts) prior_intercept else NULL,
                            prior_trt = prior_trt,
                            prior_het = if (trt_effects == "random") prior_het else NULL,
                            prior_het_type = if (trt_effects == "random") prior_het_type else NULL,
                            prior_reg = if (!is.null(regression) && !is_only_offset(regression)) prior_reg else NULL,
                            prior_aux = if (has_aux) prior_aux else NULL))

  if (inherits(network, "mlnmr_data")) class(out) <- c("stan_mlnmr", "stan_nma")
  else class(out) <- "stan_nma"

  return(out)
}


#' @param ipd_x Design matrix for IPD studies
#' @param ipd_y Outcome data frame for IPD studies
#' @param agd_arm_x  Design matrix for AgD studies (arm-based)
#' @param agd_arm_y  Outcome data frame for AgD studies (arm-based)
#' @param agd_contrast_x  Design matrix for AgD studies (contrast-based)
#' @param agd_contrast_y  Outcome data frame for AgD studies (contrast-based)
#' @param agd_contrast_Sigma List of covariance matrices for contrast-based data
#' @param n_int Number of numerical integration points used
#' @param ipd_offset Vector of offset values for IPD
#' @param agd_arm_offset Vector of offset values for AgD (arm-based)
#' @param agd_contrast_offset Vector of offset values for AgD (contrast-based)
#' @param RE_cor Random effects correlation matrix, when `trt_effects = "random"`
#' @param which_RE Random effects design vector, when `trt_effects = "random"`
#'
#' @noRd
nma.fit <- function(ipd_x, ipd_y,
                    agd_arm_x, agd_arm_y,
                    agd_contrast_x, agd_contrast_y, agd_contrast_Sigma,
                    n_int,
                    ipd_offset = NULL, agd_arm_offset = NULL, agd_contrast_offset = NULL,
                    trt_effects = c("fixed", "random"),
                    RE_cor = NULL,
                    which_RE = NULL,
                    likelihood = NULL,
                    link = NULL,
                    ...,
                    prior_intercept,
                    prior_trt,
                    prior_het,
                    prior_het_type = c("sd", "var", "prec"),
                    prior_reg,
                    prior_aux,
                    QR = FALSE,
                    adapt_delta = NULL,
                    int_thin = 100L) {

  if (missing(ipd_x)) ipd_x <- NULL
  if (missing(ipd_y)) ipd_y <- NULL
  if (missing(agd_arm_x)) agd_arm_x <- NULL
  if (missing(agd_arm_y)) agd_arm_y <- NULL
  if (missing(agd_contrast_x)) agd_contrast_x <- NULL
  if (missing(agd_contrast_y)) agd_contrast_y <- NULL
  if (missing(agd_contrast_Sigma)) agd_contrast_Sigma <- NULL

  # Check available x and y
  if (xor(is.null(ipd_x), is.null(ipd_y)))
    abort("`ipd_x` and `ipd_y` should both be present or both NULL.")
  if (xor(is.null(agd_arm_x), is.null(agd_arm_y)))
    abort("`agd_arm_x` and `agd_arm_y` should both be present or both NULL.")
  if (xor(is.null(agd_contrast_x), is.null(agd_contrast_y)) || xor(is.null(agd_contrast_x), is.null(agd_contrast_Sigma)))
    abort("`agd_contrast_x`, `agd_contrast_y`, `agd_contrast_Sigma` should all be present or all NULL.")

  has_ipd <- !is.null(ipd_x) && !is.null(ipd_y)
  has_agd_arm <- !is.null(agd_arm_x) && !is.null(agd_arm_y)
  has_agd_contrast <-  !is.null(agd_contrast_x) && !is.null(agd_contrast_y) && !is.null(agd_contrast_Sigma)

  # Ignore n_int if no AgD
  if (!has_agd_arm && !has_agd_contrast) n_int <- 1
  # Check n_int
  if (!rlang::is_scalar_integerish(n_int) ||
      n_int < 1) abort("`n_int` should be an integer >= 1.")

  # Check design matrices, outcomes
  if (has_ipd) {
    if (!is.matrix(ipd_x) || !is.numeric(ipd_x))
      abort("`ipd_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(ipd_y, is.numeric)))
      abort("`ipd_y` should be numeric outcome data.")
    if (nrow(ipd_x) != nrow(ipd_y))
      abort("Number of rows in `ipd_x` and `ipd_y` do not match.")
  }
  if (has_agd_arm) {
    if (!is.matrix(agd_arm_x) || !is.numeric(agd_arm_x))
      abort("`agd_arm_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_arm_y, is.numeric)))
      abort("`agd_arm_y` should be numeric outcome data.")
    if (nrow(agd_arm_x) != nrow(agd_arm_y) * n_int)
      abort("Number of rows in `agd_arm_x`, `agd_arm_y`, and `n_int` do not match.")
  }
  if (has_agd_contrast) {
    if (!is.matrix(agd_contrast_x) || !is.numeric(agd_contrast_x))
      abort("`agd_contrast_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_contrast_y, is.numeric)))
      abort("`agd_contrast_y` should be numeric outcome data.")
    if (nrow(agd_contrast_x) != nrow(agd_contrast_y) * n_int)
      abort("Number of rows in `agd_contrast_x`, `agd_contrast_y`, and `n_int` do not match.")
    if (!is.list(agd_contrast_Sigma) || any(purrr::map_lgl(agd_contrast_Sigma, ~!is.numeric(.))))
      abort("`agd_contrast_Sigma` should be a list of covariance matrices, of length equal to the number of AgD (contrast-based) studies.")
  }

  # Check offsets
  has_ipd_offset <- !is.null(ipd_offset)
  has_agd_arm_offset <- !is.null(agd_arm_offset)
  has_agd_contrast_offset <- !is.null(agd_contrast_offset)

  if (has_ipd_offset && !has_ipd)
    abort("`ipd_offset` given but not `ipd_x` and `ipd_y`.")
  if (has_agd_arm_offset && !has_agd_arm)
    abort("`agd_arm_offset` given but not `agd_arm_x` and `agd_arm_y`.")
  if (has_agd_contrast_offset && !has_agd_contrast)
    abort("`agd_contrast_offset` given but not `agd_contrast_x` and `agd_contrast_y`.")

  has_offsets <- any(has_ipd_offset, has_agd_arm_offset, has_agd_contrast_offset)
  if (has_offsets &&
      !all(has_ipd_offset[has_ipd],
           has_agd_arm_offset[has_agd_arm],
           has_agd_contrast_offset[has_agd_contrast]))
    abort("Offsets provided for some data sources but not all.")

  if (has_ipd_offset && (!is.numeric(ipd_offset) || length(ipd_offset) != nrow(ipd_x)))
    abort("`ipd_offset` should be a numeric vector with length matching `ipd_x`, `ipd_y`")
  if (has_agd_arm_offset && (!is.numeric(agd_arm_offset) || length(agd_arm_offset) != nrow(agd_arm_x)))
    abort("`agd_arm_offset` should be a numeric vector with length matching `agd_arm_x`, `agd_arm_y`")
  if (has_agd_contrast_offset && (!is.numeric(agd_contrast_offset) || length(agd_contrast_offset) != nrow(agd_contrast_x)))
    abort("`agd_contrast_offset` should be a numeric vector with length matching `agd_contrast_x`, `agd_contrast_y`")

  # Check matching X column names and dimensions if more than one present
  if ((has_ipd &&
       ((has_agd_arm && !identical(colnames(ipd_x), colnames(agd_arm_x))) ||
        (has_agd_contrast && !identical(colnames(ipd_x), colnames(agd_contrast_x))))) ||
      (has_agd_arm &&
       (has_agd_contrast && !identical(colnames(agd_arm_x), colnames(agd_contrast_x)))))
    abort("Non-matching columns in *_x matrices.")

  # Check model arguments
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")

  likelihood <- check_likelihood(likelihood)
  link <- check_link(link, likelihood)

  # When are priors on auxiliary parameters required?
  has_aux <- (likelihood == "normal" && has_ipd) || likelihood == "ordered"

  # Check priors
  if (!inherits(prior_intercept, "nma_prior")) abort("`prior_intercept` should be a prior distribution, see ?priors.")
  if (!inherits(prior_trt, "nma_prior")) abort("`prior_trt` should be a prior distribution, see ?priors.")
  if (trt_effects == "random" && !inherits(prior_het, "nma_prior")) abort("`prior_het` should be a prior distribution, see ?priors.")
  if (!inherits(prior_reg, "nma_prior")) abort("`prior_reg` should be a prior distribution, see ?priors.")
  if (has_aux && !inherits(prior_aux, "nma_prior")) abort("`prior_aux` should be a prior distribution, see ?priors.")

  prior_het_type <- rlang::arg_match(prior_het_type)

  # Dummy RE prior for FE model, not used but requested by Stan data
  if (trt_effects == "fixed") prior_het <- half_normal(1)

  # Check other args
  if (!rlang::is_bool(QR)) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_scalar_integerish(int_thin) ||
      int_thin < 1) abort("`int_thin` should be an integer >= 1.")

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!is.numeric(adapt_delta) ||
             length(adapt_delta) > 1 ||
             adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Pull study and treatment details from *_x
  if (has_ipd) x_names <- colnames(ipd_x)
  else if (has_agd_arm) x_names <- colnames(agd_arm_x)
  else if (has_agd_contrast) x_names <- colnames(agd_contrast_x)

  col_study <- grepl("^\\.study[^:]+$", x_names)
  col_trt <- grepl("^(\\.trt|\\.contr)[^:]+$", x_names)
  col_reg <- !col_study & !col_trt

  n_trt <- sum(col_trt) + 1

  get_study <- function(x) which(x == 1)
  get_trt <- function(x, v = 1) if (any(x == v)) which(x == v) + 1 else 1

  if (has_ipd) {
    ipd_s_t_all <- dplyr::tibble(.study = unname(apply(ipd_x[, col_study, drop = FALSE], 1, get_study)),
                                 .trt = unname(apply(ipd_x[, col_trt, drop = FALSE], 1, get_trt)))
    ipd_s_t <- dplyr::distinct(ipd_s_t_all) %>% dplyr::mutate(.arm = 1:dplyr::n())
    ipd_arm <-  dplyr::left_join(ipd_s_t_all, ipd_s_t, by = c(".study", ".trt")) %>% dplyr::pull(.data$.arm)
    ipd_study <- ipd_s_t$.study
    ipd_trt <- ipd_s_t$.trt
    narm_ipd <- max(ipd_arm)
    ni_ipd <- nrow(ipd_x)
  } else {
    ipd_study <- ipd_trt <- ipd_arm <- numeric()
    ni_ipd <- 0
    narm_ipd <- 0
  }

  if (has_agd_arm) {
    ni_agd_arm <- nrow(agd_arm_y)
    aa1 <- 0:(ni_agd_arm - 1)*n_int + 1
    agd_arm_study <- apply(agd_arm_x[aa1, col_study, drop = FALSE], 1, get_study)
    agd_arm_trt <- apply(agd_arm_x[aa1, col_trt, drop = FALSE], 1, get_trt)
  } else {
    agd_arm_study <- agd_arm_trt <- numeric()
    ni_agd_arm <- 0
  }

  if (has_agd_contrast) {
    ni_agd_contrast <- nrow(agd_contrast_y)
    ac1 <- 0:(ni_agd_contrast - 1)*n_int + 1
    agd_contrast_trt <- apply(agd_contrast_x[ac1, col_trt, drop = FALSE], 1, get_trt)
    agd_contrast_trt_b <- apply(agd_contrast_x[ac1, col_trt, drop = FALSE], 1, get_trt, v = -1)

    # Get number of contrast-based studies from length of Sigma list
    ns_agd_contrast <- length(agd_contrast_Sigma)

    # Construct block-diagonal contrast covariance matrix
    Sigma <- as.matrix(Matrix::bdiag(agd_contrast_Sigma))

    if (nrow(Sigma) != ni_agd_contrast)
      abort("Dimensions of `agd_contrast_Sigma` covariance matrices do not match the contrast-based data.")
  } else {
    agd_contrast_trt <- agd_contrast_trt_b <- numeric()
    Sigma <- matrix(1, 1, 1)
    ni_agd_contrast <- 0
    ns_agd_contrast <- 0
  }

  # Set up random effects
  if (trt_effects == "random") {
    narm <- narm_ipd + ni_agd_arm + ni_agd_contrast
    if (!is.null(which_RE)) {
      if (!rlang::is_integerish(which_RE) ||
          any(which_RE < 0) ||
          is.matrix(which_RE))
        abort("`which_RE` should be an integer vector.")
      if (length(which_RE) != narm)
        abort(glue::glue("Length of `which_RE` does not match the number of arms/contrasts.\n",
                         "Expecting length {narm}, instead length {length(which_RE)}."))
    } else {
      abort("Specify `which_RE` when trt_effects = 'random'.")
    }

    nRE <- sum(which_RE != 0)

    if (!is.null(RE_cor)) {
      if (!is.matrix(RE_cor) || !is.numeric(RE_cor))
        abort("`RE_cor` should be a numeric matrix.")
      if (any(dim(RE_cor) != c(nRE, nRE)))
        abort(glue::glue("Dimensions of `RE_cor` do not match the number of random effects.\n",
                         "Expecting [{nRE} x {nRE}], instead [{nrow(RE_cor)} x {ncol(RE_cor)}]."))
    } else {
      abort("Specify `RE_cor` when trt_effects = 'random'.")
    }
  } else {
    RE_cor <- matrix(1, 1, 1)
    which_RE <- integer()
  }

  # Make full design matrix
  X_all <- rbind(ipd_x, agd_arm_x, agd_contrast_x)

  # Make sure columns of X_all are in correct order (study, trt, regression terms)
  X_all <- cbind(X_all[, col_study, drop = FALSE],
                 X_all[, col_trt, drop = FALSE],
                 X_all[, col_reg, drop = FALSE])

  # Take thin QR decomposition if QR = TRUE
  if (QR) {
    X_all_qr <- qr(X_all)
    X_all_Q <- qr.Q(X_all_qr) * sqrt(nrow(X_all) - 1)
    X_all_R <- qr.R(X_all_qr)[, sort.list(X_all_qr$pivot)] / sqrt(nrow(X_all) - 1)
    X_all_R_inv <- solve(X_all_R)
  }

  # Set common Stan data
  standat <- list(
    # Constants
    ns_ipd = length(unique(ipd_study)),
    ni_ipd = ni_ipd,
    ns_agd_arm = length(unique(agd_arm_study)),
    ni_agd_arm = ni_agd_arm,
    ns_agd_contrast = ns_agd_contrast,
    ni_agd_contrast = ni_agd_contrast,
    nt = n_trt,
    nint = n_int,
    nX = ncol(X_all),
    int_thin = int_thin,
    # Study and treatment details
    narm_ipd = narm_ipd,
    ipd_arm = ipd_arm,
    ipd_trt = ipd_trt,
    agd_arm_trt = agd_arm_trt,
    agd_contrast_trt = as.array(agd_contrast_trt),
    agd_contrast_trt_b = as.array(agd_contrast_trt_b),
    agd_contrast_y = if (has_agd_contrast) as.array(agd_contrast_y$.y) else numeric(),
    agd_contrast_Sigma = Sigma,
    # ipd_study = ipd_study,
    # agd_arm_study = agd_arm_study,
    # agd_contrast_study = agd_contrast_study,
    # Random effects
    RE = switch(trt_effects, fixed = 0, random = 1),
    RE_cor = RE_cor,
    which_RE = which_RE,
    # Design matrix or QR decomposition
    QR = QR,
    X = if (QR) X_all_Q else X_all,
    R_inv = if (QR) X_all_R_inv else matrix(0, 0, 0),
    # Offsets
    has_offset = has_offsets,
    offsets = if (has_offsets) as.array(c(ipd_offset, agd_arm_offset, agd_contrast_offset)) else numeric()
    )

  # Add priors
  standat <- purrr::list_modify(standat,
    !!! prior_standat(prior_intercept, "prior_intercept",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_trt, "prior_trt",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_reg, "prior_reg",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_het, "prior_het",
                      valid = c("Normal", "half-Normal", "log-Normal",
                                "Cauchy",  "half-Cauchy",
                                "Student t", "half-Student t",
                                "Exponential", "flat (implicit)")),
    prior_het_type = switch(prior_het_type,
                            sd = 1, var = 2, prec = 3)
    )

  # Standard pars to monitor
  pars <- c("mu", "beta", "d",
            "log_lik", "resdev",
            "lp__")

  if (has_ipd) pars <- c(pars, "fitted_ipd")
  if (has_agd_arm) pars <- c(pars, "fitted_agd_arm")
  if (has_agd_contrast) pars <- c(pars, "fitted_agd_contrast")

  # Monitor heterogeneity SD and study deltas if RE model
  if (trt_effects == "random") {
    pars <- c(pars, "tau", "delta")
  }
  # Monitor cumulative integration error if using numerical integration
  if (n_int > 1) {
    if (has_agd_arm) pars <- c(pars, "theta_bar_cum_agd_arm")
    if (has_agd_contrast) pars <- c(pars, "theta_bar_cum_agd_contrast")
  }

  # Set adapt_delta, but respect other control arguments if passed in ...
  stanargs <- list(...)
  if ("control" %in% names(stanargs))
    stanargs$control <- purrr::list_modify(stanargs$control, adapt_delta = adapt_delta)
  else
    stanargs$control <- list(adapt_delta = adapt_delta)

  # Call Stan model for given likelihood

  # -- Normal likelihood
  if (likelihood == "normal") {

    # Dummy prior for IPD variance when no IPD - not used, but requested by Stan data
    if (!has_ipd) prior_aux <- half_normal(1)

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_y = if (has_ipd) ipd_y$.y else numeric(),
      agd_arm_y = if (has_agd_arm) agd_arm_y$.y else numeric(),
      agd_arm_se = if (has_agd_arm) agd_arm_y$.se else numeric(),

      # Add prior for auxiliary parameter - individual-level variance
      !!! prior_standat(prior_aux, "prior_aux",
                        valid = c("Normal", "half-Normal", "log-Normal",
                                  "Cauchy",  "half-Cauchy",
                                  "Student t", "half-Student t",
                                  "Exponential", "flat (implicit)")),

      # Specify link
      link = switch(link, identity = 1, log = 2)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$normal,
                                   data = standat,
                                   pars = c(pars, "sigma"))

  # -- Bernoulli/binomial likelihood (one parameter)
  } else if (likelihood %in% c("bernoulli", "binomial")) {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$binomial_1par,
                                   data = standat,
                                   pars = pars)

  # -- Bernoulli/binomial likelihood (two parameter)
  } else if (likelihood %in% c("bernoulli2", "binomial2")) {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$binomial_2par,
                                   data = standat,
                                   pars = c(pars, "theta2_bar_cum"))

  # -- Poisson likelihood
  } else if (likelihood == "poisson") {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      ipd_E = if (has_ipd) ipd_y$.E else numeric(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_E = if (has_agd_arm) agd_arm_y$.E else numeric(),

      # Specify link
      link = switch(link, log = 1)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$poisson,
                                   data = standat,
                                   pars = pars)

  # -- Ordered multinomial likelihood
  } else if (likelihood == "ordered") {

    if (has_ipd) {
      # Determine number of categories
      ncat <- ncol(ipd_y$.r)
      # Stan model takes IPD as an integer vector of category numbers
      ipd_r_int <- apply(ipd_y$.r == 1, 1, which)
      # Determine which categories are present
      ipd_cat <- t(apply(ipd_y$.r, 1,
                         function(x) {
                           cs <- which(!is.na(x))
                           c(cs, rep(0, ncat - length(cs)))
                         }))
      ipd_ncat <- rowSums(ipd_cat > 0)
    }

    if (has_agd_arm) {
      # Determine number of categories
      if (!has_ipd) ncat <- ncol(agd_arm_y$.r)
      # Determine which categories are present
      agd_arm_cat <- t(apply(agd_arm_y$.r, 1,
                         function(x) {
                           cs <- which(!is.na(x))
                           c(cs, rep(0, ncat - length(cs)))
                         }))
      agd_arm_ncat <- rowSums(agd_arm_cat > 0)
      # Replace missing category counts with 0 (these will drop out of the likelihood)
      agd_arm_r <- tidyr::replace_na(agd_arm_y$.r, 0)
      agd_arm_n <- rowSums(agd_arm_y$.r, na.rm = TRUE)
    }

    if (!has_ipd && !has_agd_arm) {
      abort("No IPD or AgD (arm-based) in the network. Cannot fit ordered model to contrast data only.")
    }

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ncat = ncat,

      ipd_r = if (has_ipd) ipd_r_int else integer(),
      ipd_cat = if (has_ipd) ipd_cat else matrix(0, 0, ncat),
      ipd_ncat = if (has_ipd) ipd_ncat else integer(),

      agd_arm_r = if (has_agd_arm) agd_arm_r else matrix(0, 0, ncat),
      agd_arm_n = if (has_agd_arm) agd_arm_n else integer(),
      agd_arm_cat = if (has_agd_arm) agd_arm_cat else matrix(0, 0, ncat),
      agd_arm_ncat = if (has_agd_arm) agd_arm_ncat else integer(),

      # Add prior for auxiliary parameters - latent cutoffs
      !!! prior_standat(prior_aux, "prior_aux",
                        valid = c("Normal", "half-Normal", "log-Normal",
                                  "Cauchy",  "half-Cauchy",
                                  "Student t", "half-Student t",
                                  "Exponential", "flat (implicit)")),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$ordered_multinomial,
                                   data = standat,
                                   pars = c(pars, "cc"))

  } else {
    abort(glue::glue('"{likelihood}" likelihood not supported.'))
  }

  stanfit <- do.call(rstan::sampling, stanargs)

  # Set readable parameter names in the stanfit object
  fnames_oi <- stanfit@sim$fnames_oi
  x_names_sub <- gsub("^(\\.study|\\.trt|\\.trtclass|\\.contr)", "", x_names)

  fnames_oi[grepl("^mu\\[[0-9]+\\]$", fnames_oi)] <- paste0("mu[", x_names_sub[col_study], "]")
  fnames_oi[grepl("^d\\[[0-9]+\\]$", fnames_oi)] <- paste0("d[", x_names_sub[col_trt], "]")
  fnames_oi[grepl("^beta\\[[0-9]+\\]$", fnames_oi)] <- paste0("beta[", x_names[col_reg], "]")
  fnames_oi <- gsub("tau[1]", "tau", fnames_oi, fixed = TRUE)

  if (likelihood == "ordered") {
    if (has_ipd) l_cat <- colnames(ipd_y$.r)[-1]
    else if (has_agd_arm) l_cat <- colnames(agd_arm_y$.r)[-1]
    fnames_oi[grepl("^cc\\[[0-9]+\\]$", fnames_oi)] <- paste0("cc[", l_cat, "]")
  }

  stanfit@sim$fnames_oi <- fnames_oi

  return(stanfit)
}

#' Random effects structure
#'
#' Use `RE_cor` to generate the random effects correlation matrix, under the
#' assumption of common heterogeneity variance (i.e. all within-study
#' correlations are 0.5). Use `which_RE` to return a vector of IDs for the RE
#' deltas (0 means no RE delta on this arm).
#'
#' @param study A vector of study IDs (integer, character, or factor)
#' @param trt A factor vector of treatment codes (or coercible as such), with
#'   first level indicating the reference treatment
#' @param contrast A logical vector, of the same length as `study` and `trt`,
#'   indicating whether the corresponding data are in contrast rather than arm
#'   format.
#' @param type Character string, whether to generate RE structure under the
#'   "reference treatment" parameterisation, or the "baseline shift"
#'   parameterisation.
#'
#' @return For `RE_cor()`, a correlation matrix of dimension equal to the number
#'   of random effects deltas (excluding those that are set equal to zero).
#'
#'   For `which_RE()`, an integer vector of IDs indexing the rows and columns of
#'   the correlation matrix returned by `RE_cor()`.
#' @export
#' @aliases RE_cor
#' @rdname random_effects
#'
#' @examples
#' RE_cor(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)))
#' RE_cor(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)), type = "blshift")
RE_cor <- function(study, trt, contrast, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (!is.logical(contrast) || is.matrix(contrast))
    abort("`contrast` must be a logical vector.")
  if (length(study) != length(trt) || length(study) != length(contrast))
    abort("`study`, `trt`, and `contrast` must be the same length.")
  type <- rlang::arg_match(type)


  reftrt <- levels(trt)[1]
  if (type == "reftrt") {
    # Treat contrast rows as non ref trt arms (since they always have REs)
    nonref <- trt != reftrt | contrast
    nRE <- sum(nonref)  # RE for each non ref trt arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[nonref]
    trt <- trt[nonref]
  } else if (type == "blshift") {
    # Consider baseline arm to be first by treatment order
    # Treat contrast rows as non baseline arms (since they always have REs)
    nonbl <- tibble::tibble(study, trt, contrast) %>%
      dplyr::group_by(.data$study) %>%
      dplyr::mutate(nonbl = .data$contrast | .data$trt != sort(.data$trt)[1] | (duplicated(.data$trt) & .data$trt != reftrt)) %>%
      dplyr::pull(.data$nonbl)
    nRE <- sum(nonbl)  # RE for each non baseline arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[nonbl]
    trt <- trt[nonbl]
  }

  diag(Rho) <- 1
  for (i in 1:(nRE - 1)) {
    for (j in (i + 1):nRE) {
      if (study[i] == study[j]) {
        Rho[i, j] <- 0.5
        Rho[j, i] <- 0.5
      }
    }
  }

  return(Rho)
}

#' @rdname random_effects
#' @aliases which_RE
#' @export
#' @examples
#' which_RE(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)))
#' which_RE(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)), type = "blshift")
which_RE <- function(study, trt, contrast, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (!is.logical(contrast) || is.matrix(contrast))
    abort("`contrast` must be a logical vector.")
  if (length(study) != length(trt) || length(study) != length(contrast))
    abort("`study`, `trt`, and `contrast` must be the same length.")
  type <- rlang::arg_match(type)

  n_i <- length(study)
  id <- rep(0, n_i)

  reftrt <- levels(trt)[1]
  if (type == "reftrt") {
    trt_nonref <- trt != reftrt | contrast
    id[trt_nonref] <- 1:sum(trt_nonref)
  } else if (type == "blshift") {
    # Consider baseline arm to be first by treatment order
    non_bl <- tibble::tibble(study, trt, contrast) %>%
      dplyr::group_by(.data$study) %>%
      dplyr::mutate(non_bl = .data$contrast | .data$trt != sort(.data$trt)[1] | (duplicated(.data$trt) & .data$trt != reftrt)) %>%
      dplyr::pull(.data$non_bl)
    id[non_bl] <- 1:sum(non_bl)
  }

  return(id)
}

#' Check likelihood function, or provide default value
#'
#' @param x likelihood type as string
#' @param outcome outcome types as named list (ipd, agd_arm, agd_contrast)
#'
#' @noRd
check_likelihood <- function(x, outcome) {
  valid_lhood <- list(binary = c("bernoulli", "bernoulli2"),
                      count = c("binomial", "binomial2"),
                      rate = "poisson",
                      continuous = "normal",
                      ordered = "ordered")

  if (missing(outcome)) valid_lhood <- unlist(valid_lhood)
  else if (!is.na(outcome$ipd)) {
    valid_lhood <- valid_lhood[[outcome$ipd]]
    otype <- outcome$ipd
  } else if (!is.na(outcome$agd_arm)) {
    valid_lhood <- valid_lhood[[outcome$agd_arm]]
    otype <- outcome$agd_arm
  } else if (!is.na(outcome$agd_contrast)) {
    valid_lhood <- valid_lhood[[outcome$agd_contrast]]
    otype <- outcome$agd_contrast
  }
  else abort("No outcomes specified in network data.")

  # Choose default likelihood if NULL
  if (is.null(x) && missing(outcome)) abort("Please specify a suitable likelihood.")
  else if (is.null(x)) {
    x <- valid_lhood[1]

  # Check valid option if given
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_lhood) {
    abort(glue::glue("`likelihood` should be a character string specifying a valid likelihood.\n",
                     "Suitable options for {otype} outcomes are currently: ",
                     glue::glue_collapse(glue::double_quote(valid_lhood), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

#' Check link function, or provide default value
#'
#' @param x link function, as string
#' @param lik likelihood, as string
#'
#' @noRd
check_link <- function(x, lik) {
  valid_link <- list(normal = c("identity", "log"),
                     bernoulli = c("logit", "probit", "cloglog"),
                     bernoulli2 = c("logit", "probit", "cloglog"),
                     binomial = c("logit", "probit", "cloglog"),
                     binomial2 = c("logit", "probit", "cloglog"),
                     poisson = "log",
                     ordered = c("logit", "probit", "cloglog"))[[lik]]

  if (is.null(x)) {
    x <- valid_link[1]
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_link) {
    abort(glue::glue("`link` should be a character string specifying a valid link function.\n",
                     "Suitable options for a {lik} likelihood are currently: ",
                     glue::glue_collapse(glue::double_quote(valid_link), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

#' Inverse link functions
#'
#' @param x Linear predictor values
#' @param link Character string specifying link function
#' @param ... Other parameters passed to link function
#'
#' @noRd
inverse_link <- function(x, link = c("identity", "log", "logit", "probit", "cloglog"), ...) {
  link <- rlang::arg_match(link)

  out <-
    if (link == "identity") x
    else if (link == "log") exp(x)
    else if (link == "logit") plogis(x, ...)
    else if (link == "probit") pnorm(x, ...)
    else if (link == "cloglog") 1 - exp(-exp(x))

  return(out)
}

#' Get scale of outcome / linear predictor for reporting and plotting
#'
#' @param likelihood String, giving likelihood
#' @param link String, giving link function
#' @param measure String, specifying whether relative or absolute scale required
#' @param type String, specifying whether link scale or response scale required
#'
#' @return String giving the scale name, e.g. "log Odds Ratio"
#' @noRd
get_scale_name <- function(likelihood = c("normal", "bernoulli", "bernoulli2",
                                          "binomial", "binomial2", "poisson",
                                          "ordered"),
                           link = c("identity", "log", "logit", "probit", "cloglog"),
                           measure = c("relative", "absolute"),
                           type = c("link", "response")) {

  likelihood <- rlang::arg_match(likelihood)
  link <- rlang::arg_match(link)
  measure <- rlang::arg_match(measure)
  type <- rlang::arg_match(type)

  check_link(link, likelihood)

  if (likelihood == "normal") {

    if (link == "identity") {
      if (measure == "relative") {
        out <- "Relative Effect"
      } else if (measure == "absolute") {
        out <- "Absolute Effect"
      }
    } else if (link == "log") {
      if (measure == "relative") {
        if (type == "link") out <- "Relative Effect (on log scale)"
        else out <- "Relative Effect (on response scale)"
      } else if (measure == "absolute") {
        if (type == "link") out <- "Absolute Effect (on log scale)"
        else out <- "Absolute Effect (on response scale)"
      }
    }

  } else if (likelihood %in% c("bernoulli", "bernoulli2", "binomial", "binomial2", "ordered")) {

    if (link == "logit") {
      if (measure == "relative") {
        if (type == "link") out <- "log Odds Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "log Odds"
        else out <- "Probability"
      }
    } else if (link == "probit") {
      if (measure == "relative") {
        if (type == "link") out <- "Probit Difference"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "Probit Probability"
        else out <- "Probability"
      }
    } else if (link == "cloglog") {
      if (measure == "relative") {
        if (type == "link") out <- "log Hazard Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "cloglog Probability"
        else out <- "Probability"
      }
    }

  } else if (likelihood == "poisson") {

    if (link == "log") {
      if (measure == "relative") {
        if (type == "link") out <- "log Rate Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "log Rate"
        else out <- "Rate"
      }
    }

  } else {
    out <- ""
  }
  return(out)
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
    continuous = c(".y", ".se"),
    ordered = ".r"
  )[[o_type]]

  return(
    dplyr::select(x, dplyr::one_of(intersect(o_vars, colnames(x))))
  )
}

#' Construct NMA formula
#'
#' @param regression Regression formula, or NULL
#' @param consistency Consistency/inconsistency model (character string)
#' @param classes Classes present? TRUE / FALSE
#' @param class_interactions Class interaction specification (character string)
#'
#' @return A formula
#' @noRd
make_nma_formula <- function(regression,
                             consistency = c("consistency", "nodesplit", "ume"),
                             classes,
                             class_interactions = c("common", "exchangeable", "independent")
                             ) {

  if (!is.null(regression) && !rlang::is_formula(regression)) abort("`regression` is not a formula")
  consistency <- rlang::arg_match(consistency)
  if (!rlang::is_bool(classes)) abort("`classes` should be TRUE or FALSE")
  if (classes && !is.null(regression)) class_interactions <- rlang::arg_match(class_interactions)

  if (!is.null(regression)) {

    # Set up treatment classes
    if (classes) {

      if (class_interactions == "common") {
        nma_formula <- do.call("substitute",
                               list(regression,
                                    list(.trt = quote(.trtclass))))

        # Remove any main effect of .trtclass, e.g. if user specified var*.trt
        nma_formula <- update.formula(nma_formula, ~. - .trtclass)

      } else if (class_interactions == "exchangeable") {
        abort('Exchangeable treatment class interactions (class_interactions = "exchangeable") not yet supported.')
      } else {
        nma_formula <- regression
      }
    } else {
      nma_formula <- regression
    }

    if (consistency == "ume") {
      nma_formula <- update.formula(nma_formula, ~-1 + .study + .contr + .)
    } else {
      nma_formula <- update.formula(nma_formula, ~-1 + .study + .trt + .)
    }
  } else {
    if (consistency == "ume") {
      nma_formula <- ~-1 + .study + .contr
    } else {
      nma_formula <- ~-1 + .study + .trt
    }
  }
}

#' Construct NMA design matrix
#'
#' @param nma_formula NMA formula, returned by [make_nma_formula()]
#' @param ipd,agd_arm,agd_contrast Data frames
#' @param agd_contrast_bl Logical vector identifying baseline rows for contrast
#'   data
#' @param xbar Named numeric vector of centering values, or NULL
#' @param consistency Consistency/inconsistency model (character string)
#' @param classes Classes present? TRUE / FALSE
#' @param class_interactions Class interaction specification (character string)
#'
#' @return A named list of three matrices: X_ipd, X_agd_arm, X_agd_contrast; and
#'   three vectors of offsets: offset_ipd, offset_agd_arm, offset_agd_contrast.
#' @noRd
make_nma_model_matrix <- function(nma_formula,
                                  dat_ipd = tibble::tibble(),
                                  dat_agd_arm = tibble::tibble(),
                                  dat_agd_contrast = tibble::tibble(),
                                  agd_contrast_bl = logical(),
                                  xbar = NULL,
                                  consistency = c("consistency", "nodesplit", "ume"),
                                  classes = FALSE) {
  # Checks
  if (!rlang::is_formula(nma_formula)) abort("`nma_formula` is not a formula")
  stopifnot(is.data.frame(dat_ipd),
            is.data.frame(dat_agd_arm),
            is.data.frame(dat_agd_contrast))
  consistency <- rlang::arg_match(consistency)
  if (!rlang::is_bool(classes)) abort("`classes` should be TRUE or FALSE")
  if (nrow(dat_agd_contrast) && !rlang::is_logical(agd_contrast_bl, n = nrow(dat_agd_contrast)))
    abort("`agd_contrast_bl` should be a logical vector of length nrow(agd_contrast)")
  if (!is.null(xbar) && (
        !(rlang::is_double(xbar) || rlang::is_integer(xbar)) || !rlang::is_named(xbar)))
    abort("`xbar` should be a named numeric vector")

  if (!consistency %in% c("consistency", "ume")) {
    abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
  }

  .has_ipd <- if (nrow(dat_ipd)) TRUE else FALSE
  .has_agd_arm <- if (nrow(dat_agd_arm)) TRUE else FALSE
  .has_agd_contrast <- if (nrow(dat_agd_contrast)) TRUE else FALSE

  # Sanitise factors
  if (.has_ipd) {
    dat_ipd <- dplyr::mutate_at(dat_ipd,
      .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
      .funs = fct_sanitise)
  }
  if (.has_agd_arm) {
    dat_agd_arm <- dplyr::mutate_at(dat_agd_arm,
                            .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
                            .funs = fct_sanitise)
  }
  if (.has_agd_contrast) {
    dat_agd_contrast <- dplyr::mutate_at(dat_agd_contrast,
                            .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
                            .funs = fct_sanitise)

    # Split contrast-based data into baseline and non-baseline arms
    dat_agd_contrast_bl <- dat_agd_contrast[agd_contrast_bl, ]
    dat_agd_contrast_nonbl <- dat_agd_contrast[!agd_contrast_bl, ]
  } else {
    dat_agd_contrast_bl <- dat_agd_contrast_nonbl <- tibble::tibble()
  }

  # Define contrasts for UME model
  if (consistency == "ume") {
    # For IPD and AgD (arm-based), take the first-ordered arm as baseline
    # So the contrast sign will always be positive
    if (.has_ipd || .has_agd_arm) {
      contrs_arm <- dplyr::bind_rows(dat_ipd, dat_agd_arm) %>%
        dplyr::distinct(.data$.study, .data$.trt) %>%
        dplyr::arrange(.data$.study, .data$.trt) %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::mutate(.trt_b = dplyr::first(.data$.trt)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.contr = dplyr::if_else(.data$.trt == .data$.trt_b,
                                              "..ref..",
                                              paste0(.data$.trt, " vs. ", .data$.trt_b)),
                      .contr_sign = 1)
    } else {
      contrs_arm <- tibble::tibble()
    }

    # For AgD (contrast-based), take the specified baseline arm (with .y = NA)
    # Need to make sure to construct contrast the correct way around (d_12 instead of d_21)
    # and then make a note to change the sign if necessary
    if (.has_agd_contrast) {
      contrs_contr <- dat_agd_contrast %>%
        dplyr::distinct(.data$.study, .data$.trt) %>%  # In case integration data passed
        dplyr::group_by(.data$.study) %>%
        dplyr::mutate(.trt_b = .data$.trt[which(is.na(.data$.y))]) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.data$.study, .data$.trt, .data$.trt_b) %>%
        dplyr::mutate(.contr_sign = dplyr::if_else(as.numeric(.data$.trt) < as.numeric(.data$.trt_b), -1, 1),
                      .contr = dplyr::if_else(.data$.trt == .data$.trt_b,
                                              "..ref..",
                                              dplyr::if_else(.data$.contr_sign == 1,
                                                             paste0(.data$.trt, " vs. ", .data$.trt_b),
                                                             paste0(.data$.trt_b, " vs. ", .data$.trt))))
    } else {
      contrs_contr <- tibble::tibble()
    }

    # Make contrast info
    contrs_all <- dplyr::bind_rows(contrs_arm, contrs_contr)

    # Vector of all K(K-1)/2 possible contrast levels
    nt <- nlevels(contrs_all$.trt)
    ctr <- which(lower.tri(diag(nt)), arr.ind = TRUE)
    trt_lev <- levels(contrs_all$.trt)
    c_lev <- paste(trt_lev[ctr[, "row"]], trt_lev[ctr[, "col"]], sep = " vs. ")

    contrs_all <- dplyr::transmute(contrs_all,
      .data$.study, .data$.trt,
      .contr = forcats::fct_drop(factor(.data$.contr, levels = c("..ref..", c_lev))),
      .data$.contr_sign)

    # Join contrast info on to study data
    if (.has_ipd)
      dat_ipd <- dplyr::left_join(dat_ipd, contrs_all, by = c(".study", ".trt"))
    if (.has_agd_arm) {
      dat_agd_arm <- dplyr::left_join(dat_agd_arm, contrs_all, by = c(".study", ".trt"))
    }
    if (.has_agd_contrast) {
      dat_agd_contrast <- dplyr::left_join(dat_agd_contrast, contrs_all, by = c(".study", ".trt"))
      dat_agd_contrast_bl <- dplyr::left_join(dat_agd_contrast_bl, contrs_all, by = c(".study", ".trt"))
      dat_agd_contrast_nonbl <- dplyr::left_join(dat_agd_contrast_nonbl, contrs_all, by = c(".study", ".trt"))
    }
  }

  # Construct design matrix all together then split out, so that same dummy
  # coding is used everywhere
  dat_all <- dplyr::bind_rows(dat_ipd, dat_agd_arm, dat_agd_contrast_nonbl)

  # Check that required variables are present in each data set, and non-missing
  check_regression_data(nma_formula,
                        dat_ipd = dat_ipd,
                        dat_agd_arm = dat_agd_arm,
                        dat_agd_contrast = dat_agd_contrast)

  # Center
  if (!is.null(xbar)) {
    dat_all[, names(xbar)] <-
      purrr::map2(dat_all[, names(xbar)], xbar, ~.x - .y)

    if (.has_agd_contrast) {
      dat_agd_contrast_bl[, names(xbar)] <-
        purrr::map2(dat_agd_contrast_bl[, names(xbar)], xbar, ~.x - .y)
    }
  }

  # Drop study to factor to 1L if only one study (avoid contrasts need 2 or
  # more levels error)
  if (dplyr::n_distinct(dat_all$.study) == 1) {

    # Save study label to restore
    single_study_label <- unique(dat_all$.study)
    dat_all$.study_temp <- dat_all$.study
    dat_all$.study <- 1L

    # Fix up model formula with an intercept
    nma_formula <- update.formula(nma_formula, ~. + 1)
  } else {
    single_study_label <- NULL
  }

  # Apply NMA formula to get design matrix
  X_all <- model.matrix(nma_formula, data = dat_all)
  offsets <- model.offset(model.frame(nma_formula, data = dat_all))
  has_offset <- !is.null(offsets)

  if (!is.null(single_study_label)) {
    # Restore single study label and .study column
    colnames(X_all) <- stringr::str_replace(colnames(X_all),
                                            "^\\.study$",
                                            paste0(".study", single_study_label))
    dat_all <- dat_all %>%
      dplyr::mutate(.study = .data$.study_temp) %>%
      dplyr::select(-.data$.study_temp)

    # Drop intercept column from design matrix
    X_all <- X_all[, -1, drop = FALSE]
  }

  # Remove columns for reference level of .trtclass
  if (classes) {
    ref_class <- levels(dat_all$.trtclass)[1]
    col_trtclass_ref <- grepl(paste0(".trtclass", ref_class),
                              colnames(X_all), fixed = TRUE)
    X_all <- X_all[, !col_trtclass_ref, drop = FALSE]
  }

  # Remove columns for interactions with reference level of .trt or .trtclass
  ref_trt <- levels(dat_all$.trt)[1]
  regex_int_ref <- paste0("\\:\\.trt\\Q", ref_trt, "\\E$|^\\.trt\\Q", ref_trt, "\\E\\:")
  if (classes)
    regex_int_ref <- paste0(regex_int_ref, "|",
                            "\\:\\.trtclass\\Q", ref_class, "\\E$|^\\.trtclass\\Q", ref_class, "\\E\\:")
  col_int_ref <- grepl(regex_int_ref, colnames(X_all), perl = TRUE)
  X_all <- X_all[, !col_int_ref, drop = FALSE]

  if (consistency == "ume") {
    # Set relevant entries to +/- 1 for direction of contrast, using .contr_sign
    contr_cols <- grepl("^\\.contr", colnames(X_all))
    X_all[, contr_cols] <- sweep(X_all[, contr_cols, drop = FALSE], MARGIN = 1,
                                 STATS = dat_all$.contr_sign, FUN = "*")
  }

  if (.has_ipd) {
    X_ipd <- X_all[1:nrow(dat_ipd), , drop = FALSE]
    offset_ipd <- if (has_offset) offsets[1:nrow(dat_ipd)] else NULL
  } else {
    X_ipd <- offset_ipd <- NULL
  }

  if (.has_agd_arm) {
    X_agd_arm <- X_all[nrow(dat_ipd) + 1:nrow(dat_agd_arm), , drop = FALSE]
    offset_agd_arm <- if (has_offset) offsets[nrow(dat_ipd) + 1:nrow(dat_agd_arm)] else NULL
  } else {
    X_agd_arm <- offset_agd_arm <- NULL
  }

  if (.has_agd_contrast) {
    X_agd_contrast <- X_all[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast_nonbl), , drop = FALSE]
    offset_agd_contrast <-
      if (has_offset) {
        offsets[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast_nonbl)]
      } else {
        NULL
      }

    # Fix up single study case
    if (!is.null(single_study_label)) {
      dat_agd_contrast_bl$.study_temp <- dat_agd_contrast_bl$.study
      dat_agd_contrast_bl$.study <- 1L
    }

    # Difference out the baseline arms
    X_bl <- model.matrix(nma_formula, data = dat_agd_contrast_bl)
    if (has_offset) offset_bl <- model.offset(model.frame(nma_formula, data = dat_agd_contrast_bl))

    if (!is.null(single_study_label)) {
      # Restore single study label and .study column
      colnames(X_bl) <- stringr::str_replace(colnames(X_bl),
                                             "^\\.study$",
                                             paste0(".study", single_study_label))
      dat_agd_contrast_bl <- dat_agd_contrast_bl %>%
        dplyr::mutate(.study = .data$.study_temp) %>%
        dplyr::select(-.data$.study_temp)

      # Drop intercept column from design matrix
      X_bl <- X_bl[, -1, drop = FALSE]
    }

    # The factor levels should be the same between idat_all and
    # idat_agd_contrast_bl, so the same columns should be present in both design
    # matrices - but check anyway
    if (any(colnames(X_agd_contrast) != colnames(X_bl)))
      abort("Mismatch design matrices for baseline and non-baseline arms. Dropped factor levels?")

    if (consistency == "ume") {
      # Set relevant entries to +/- 1 for direction of contrast, using .contr_sign
      X_bl[, contr_cols] <- sweep(X_bl[, contr_cols, drop = FALSE], MARGIN = 1,
                                  STATS = dat_agd_contrast_bl$.contr_sign, FUN = "*")
    }

    # Match non-baseline rows with baseline rows by study
    bl_lookup <- vapply(dat_agd_contrast_nonbl$.study,
                        FUN = function(x) which(x == dat_agd_contrast_bl$.study),
                        FUN.VALUE = numeric(1))

    X_agd_contrast <- X_agd_contrast - X_bl[bl_lookup, , drop = FALSE]
    if (has_offset) offset_agd_contrast <- offset_agd_contrast - offset_bl[bl_lookup]

    # Remove columns for study baselines corresponding to contrast-based studies - not used
    s_contr <- unique(dat_agd_contrast$.study)
    bl_s_reg <- paste0("^\\.study(\\Q", paste0(s_contr, collapse = "\\E|\\Q"), "\\E)$")
    bl_cols <- grepl(bl_s_reg, colnames(X_agd_contrast), perl = TRUE)

    X_agd_contrast <- X_agd_contrast[, !bl_cols, drop = FALSE]
    if (.has_ipd) X_ipd <- X_ipd[, !bl_cols, drop = FALSE]
    if (.has_agd_arm) X_agd_arm <- X_agd_arm[, !bl_cols, drop = FALSE]
  } else {
    X_agd_contrast <- offset_agd_contrast <- NULL
  }

  return(list(X_ipd = X_ipd,
              X_agd_arm = X_agd_arm,
              X_agd_contrast = X_agd_contrast,
              offset_ipd = offset_ipd,
              offset_agd_arm = offset_agd_arm,
              offset_agd_contrast = offset_agd_contrast))
}

#' Extract columns used in model from data frame
#'
#' @param data Data frame
#' @param regression Regression formula or NULL
#' @param label Label for data source or NULL, used for informative errors
#'
#' @return Data frame with required columns
#' @noRd
get_model_data_columns <- function(data, regression = NULL, label = NULL) {
  if (!is.null(label)) label <- paste(" in", label)
  if (!is.null(regression)) {
    regvars <- setdiff(all.vars(regression), c(".trt", ".trtclass", ".study"))
    badvars <- setdiff(regvars, colnames(data))
    if (length(badvars)) {
      abort(
        glue::glue("Regression variable{if (length(badvars) > 1) 's' else ''} ",
                   glue::glue_collapse(glue::double_quote(badvars), sep = ", ", last = " and "),
                   " not found", label, ".")
      )
    }
    out <- dplyr::select(data, dplyr::starts_with("."), !! regvars)
  } else {
    out <- dplyr::select(data, dplyr::starts_with("."))
  }
  return(out)
}

#' Check data for regression
#'
#' Validates data for regression model. Checks that model matrix can be
#' constructed, and that there are no missing values.
#'
#' @param formula Model formula
#' @param dat_ipd,dat_agd_arm,dat_agd_contrast Data frames
#'
#' @noRd
check_regression_data <- function(formula,
                                  dat_ipd = tibble::tibble(),
                                  dat_agd_arm = tibble::tibble(),
                                  dat_agd_contrast = tibble::tibble()) {

  .has_ipd <- if (nrow(dat_ipd)) TRUE else FALSE
  .has_agd_arm <- if (nrow(dat_agd_arm)) TRUE else FALSE
  .has_agd_contrast <- if (nrow(dat_agd_contrast)) TRUE else FALSE

  # Check that required variables are present in each data set, and non-missing
  if (.has_ipd) {
    rlang::with_handlers(
      X_ipd_frame <- model.frame(formula, dat_ipd, na.action = NULL),
      error = ~abort(paste0("Failed to construct design matrix for IPD.\n", .)))

    X_ipd_has_na <- names(which(purrr::map_lgl(X_ipd_frame, ~any(is.na(.)))))
  } else {
    X_ipd_has_na <- character(0)
  }

  if (.has_agd_arm) {
    rlang::with_handlers(
      X_agd_arm_frame <- model.frame(formula, dat_agd_arm, na.action = NULL),
      error = ~abort(paste0("Failed to construct design matrix for AgD (arm-based).\n", .)))

    X_agd_arm_has_na <- names(which(purrr::map_lgl(X_agd_arm_frame, ~any(is.na(.)))))
  } else {
    X_agd_arm_has_na <- character(0)
  }

  if (.has_agd_contrast) {
    rlang::with_handlers(
      X_agd_contrast_frame <- model.frame(formula, dat_agd_contrast, na.action = NULL),
      error = ~abort(paste0("Failed to construct design matrix for AgD (contrast-based).\n", .)))

    X_agd_contrast_has_na <- names(which(purrr::map_lgl(X_agd_contrast_frame, ~any(is.na(.)))))
  } else {
    X_agd_contrast_has_na <- character(0)
  }

  dat_has_na <- c(length(X_ipd_has_na) > 0,
                  length(X_agd_arm_has_na) > 0,
                  length(X_agd_contrast_has_na) > 0)
  if (any(dat_has_na)) {
    abort(glue::glue(glue::glue_collapse(
      c("Variables with missing values in IPD: {paste(X_ipd_has_na, collapse = ', ')}.",
        "Variables with missing values in AgD (arm-based): {paste(X_agd_arm_has_na, collapse = ', ')}.",
        "Variables with missing values in AgD (contrast-based): {paste(X_agd_contrast_has_na, collapse = ', ')}."
      )[dat_has_na], sep = "\n")))
  }

  invisible()
}

#' Set prior details for Stan models
#'
#' @param x a `nma_prior` object
#' @param par character string, giving the Stan root parameter name (e.g.
#'   "prior_trt")
#' @param valid character vector, giving valid distributions
#'
#' @noRd
prior_standat <- function(x, par, valid){
  if (!inherits(x, "nma_prior")) abort("Not a `nma_prior` object.")

  dist <- x$dist

  if (!dist %in% valid)
    abort(glue::glue("Invalid `{par}`. Suitable distributions are: ",
                glue::glue_collapse(valid, sep = ", ", last = ", or ")))

  distn <- switch(dist,
                  `flat (implicit)` = 0,
                  Normal = , `half-Normal` = 1,
                  Cauchy = , `half-Cauchy` = 2,
                  `Student t` = , `half-Student t` = 3,
                  Exponential = 4,
                  `log-Normal` = 5)

  out <- purrr::list_modify(c(x), dist = distn, fun = purrr::zap())
  # Set unnecessary (NA) parameters to zero. These will be ignored by Stan, but
  # need to pass rstan checks
  out[is.na(out)] <- 0
  names(out) <- paste0(par, "_", names(out))
  return(out)
}

#' Get covariance structure contrast-based data, using se on baseline arm
#'
#' @param x A data frame of agd_contrast data
#'
#' @return A list of covariance matrices, of length equal to the number of studies
#' @noRd
make_Sigma <- function(x) {
  return(unclass(
           by(x,
              forcats::fct_inorder(forcats::fct_drop(x$.study)),
              FUN = make_Sigma_block,
              simplify = FALSE)))
}

make_Sigma_block <- function(x) {
  s_ij <- x[is.na(x$.y), ".se", drop = TRUE]^2  # Covariances
  s_ii <- x[!is.na(x$.y), ".se", drop = TRUE]^2  # Variances
  narm <- length(s_ii)
  S <- matrix(s_ij, nrow = narm, ncol = narm)
  diag(S) <- s_ii
  return(S)
}

#' Sanitise factor labels
#'
#' Remove forbidden chars from factor labels, replacing with "_"
#'
#' @param f A factor
#'
#' @noRd
fct_sanitise <- function(f) {
  forcats::fct_relabel(f, ~stringr::str_replace_all(., "[:\\[\\]]", "_"))
}

#' Check if formula only contains an offset
#'
#' @param f A one-sided formula, as used for regression argument in nma()
#' @noRd
is_only_offset <- function(f) {
  if (!rlang::is_formula(f, lhs = FALSE))
    abort("`f` should be a one-sided formula.")

  fterms <- terms(f)

  out <- !is.null(attr(fterms, "offset")) && length(attr(fterms, "factors")) == 0
  return(out)
}

#' Get EM variables from a model formula
#'
#' @param f A one-sided formula, as used for regression argument / output from make_nma_formula
#' @noRd
#' @return A character vector or character(0)
get_EM_vars <- function(f) {
  if (!rlang::is_formula(f, lhs = FALSE))
    abort("`f` should be a one-sided formula.")

  fnames <- colnames(attr(terms(f), "factor"))

  EM_regex <- "(^\\.trt(class)?\\:)|(\\:\\.trt(class)?$)"

  out <- stringr::str_subset(fnames, EM_regex)
  out <- unique(stringr::str_remove(out, EM_regex))
  return(out)
}

#' Make labels for data points
#'
#' @param study Study vector, coerced to character
#' @param trt Treatment vector, coerced to character
#' @param trt_b Baseline treatment vector if contrast-based, coerced to character
#'
#' @return Character vector of labels
#'
#' @noRd
make_data_labels <- function(study, trt, trt_b = NA) {
  tibble::tibble(study = study, trt = trt, trt_b = trt_b) %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(grp_n = dplyr::n(),
                  grp_id = 1:dplyr::n(),
                  vs_trt_b = dplyr::if_else(is.na(trt_b),
                                            NA_character_,
                                            paste0(" vs. ", trt_b)),
                  label = as.character(
                    dplyr::if_else(.data$grp_n > 1,
                                   glue::glue("{study}: {trt}{vs_trt_b}, {grp_id}", .na = ""),
                                   glue::glue("{study}: {trt}{vs_trt_b}", .na = "")))
    ) %>%
    dplyr::pull(.data$label)
}
