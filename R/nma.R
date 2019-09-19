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
#' @param center Logical scalar (default `TRUE`), whether to center the
#'   (numeric) regression terms about the overall means
#' @param agd_sample_size Column name (either bare or a string) in the AgD
#'   giving the sample size to use for calculating the global means
#' @param adapt_delta See [adapt_delta] for details
#' @param int_thin A single integer value, the thinning factor for returning
#'   cumulative estimates of integration error
#'
#' @return `nma()` returns a [stan_nma] object, `nma.fit()` returns a [stanfit]
#'   object.
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
                prior_intercept = normal(scale = 10),
                prior_trt = normal(scale = 10),
                prior_het = half_normal(scale = 5),
                prior_reg = normal(scale = 10),
                prior_aux = normal(scale = 5),
                QR = FALSE,
                center = TRUE,
                agd_sample_size,
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
  if (!is.logical(center) || length(center) > 1)
    abort("`center` should be a logical scalar (TRUE or FALSE).")
  if (!is.numeric(int_thin) ||
      length(int_thin) > 1 ||
      int_thin < 1 ||
      trunc(int_thin) != int_thin) abort("`int_thin` should be an integer >= 1.")

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!is.numeric(adapt_delta) ||
      length(adapt_delta) > 1 ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Number of numerical integration points
  # Set to 1 if no numerical integration, so that regression on summary data is possible
  n_int <- if (inherits(network, "mlnmr_data")) network$n_int else 1

  # Get design matrices and outcomes
  if (has_ipd(network)) {
    dat_ipd <- network$ipd
    y_ipd <- get_outcome_variables(dat_ipd, network$outcome$ipd)
  } else {
    dat_ipd <- tibble::tibble()
    y_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    dat_agd_arm <- network$agd_arm
    y_agd_arm <- get_outcome_variables(dat_agd_arm, network$outcome$agd_arm)

    # Set up integration variables if present
    if (inherits(network, "mlnmr_data")) {
      idat_agd_arm <- dat_agd_arm %>%
         dplyr::select(-dplyr::one_of(intersect(network$int_names, colnames(dat_agd_arm)))) %>%
         dplyr::rename_at(dplyr::vars(dplyr::starts_with(".int_")), ~gsub(".int_", "", ., fixed = TRUE)) %>%
         tidyr::unnest(!!! rlang::syms(network$int_names))
    } else {
      idat_agd_arm <- dat_agd_arm
    }
  } else {
    dat_agd_arm <- idat_agd_arm <- tibble::tibble()
    y_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    dat_agd_contrast <- network$agd_contrast
    y_agd_contrast <- get_outcome_variables(dat_agd_contrast, network$outcome$agd_contrast)

    # Set up integration variables if present
    if (inherits(network, "mlnmr_data")) {
      idat_agd_contrast <- dat_agd_contrast %>%
        dplyr::select(-dplyr::one_of(union(network$int_names, colnames(dat_agd_contrast)))) %>%
        dplyr::rename_at(dplyr::vars(dplyr::starts_with(".int_")), ~gsub(".int_", "", ., fixed = TRUE)) %>%
        tidyr::unnest(!!! rlang::syms(network$int_names))
    } else {
      idat_agd_contrast <- dat_agd_contrast
    }
  } else {
    dat_agd_contrast <- idat_agd_contrast <- tibble::tibble()
    y_agd_contrast <- NULL
  }

  # Construct design matrix all together then split out, so that same dummy
  # coding is used everywhere
  idat_all <- dplyr::bind_rows(dat_ipd, idat_agd_arm, idat_agd_contrast) %>%

    # Sanitise study and treatment factor labels (for :)
    dplyr::mutate(.study = forcats::fct_relabel(.data$.study, ~gsub(":", "_", ., fixed = TRUE)),
                  .trt = forcats::fct_relabel(.data$.trt, ~gsub(":", "_", ., fixed = TRUE)))

  # Get sample sizes for centering
  if (!is.null(regression) && center) {
    if ((has_agd_arm(network) || has_agd_contrast(network))&& missing(agd_sample_size))
      abort("Specify AgD sample size column in data `agd_sample_size` to calculate global mean for centering, or set center = FALSE.")

    if (has_agd_arm(network)) {
      N_agd_arm <- dplyr::pull(network$agd_arm, {{ agd_sample_size }})
    } else {
      N_agd_arm <- NULL
    }

  if (has_agd_contrast(network)) {

    # For contrast-based data, the "contrast" sample size is undefined -
    # expecting instead to see study sample size repeated by contrast (ditto for
    # regression terms). For centering, take the first value of N, and set
    # others to zero.
    first_then_zero <- function(x) {
      x[2:length(x)] <- 0
      return(x)
    }
      N_agd_contrast <- network$agd_contrast %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::mutate_at(dplyr::vars({{ agd_sample_size }}), first_then_zero) %>%
        dplyr::pull(network$agd_arm, {{ agd_sample_size }})
    } else {
      N_agd_contrast <- NULL
    }

  # Center numeric columns in data
    wts <- c(rep(1, nrow(dat_ipd)),
             rep(N_agd_arm / n_int, each = n_int),
             rep(N_agd_contrast / n_int, each = n_int))

    idat_all <- dplyr::mutate_if(idat_all, is.numeric, ~. - weighted.mean(., wts))
  }

  if (consistency != "consistency") {
    abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
  }

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

  # Check that required variables are present in each data set, and non-missing
  if (has_ipd(network)) {
    rlang::with_handlers(
      X_ipd_frame <- model.frame(nma_formula, dat_ipd, na.action = NULL),
      error = ~abort(paste0("Failed to construct design matrix for IPD.\n", .)))

    X_ipd_has_na <- names(which(purrr::map_lgl(X_ipd_frame, ~any(is.na(.)))))
  } else {
    X_ipd_has_na <- character(0)
  }

  if (has_agd_arm(network)) {
    rlang::with_handlers(
      X_agd_arm_frame <- model.frame(nma_formula, idat_agd_arm, na.action = NULL),
      error = ~abort(paste0("Failed to construct design matrix for AgD (arm-based).\n", .)))

    X_agd_arm_has_na <- names(which(purrr::map_lgl(X_agd_arm_frame, ~any(is.na(.)))))
  } else {
    X_agd_arm_has_na <- character(0)
  }

  if (has_agd_contrast(network)) {
    rlang::with_handlers(
      X_agd_contrast_frame <- model.frame(nma_formula, idat_agd_contrast, na.action = NULL),
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

  # Construct model matrix
  X_all <- model.matrix(nma_formula, data = idat_all)

  if (has_ipd(network)) {
    X_ipd <- X_all[1:nrow(dat_ipd), ]
  } else {
    X_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    X_agd_arm <- X_all[nrow(dat_ipd) + 1:nrow(idat_agd_arm), ]
  } else {
    X_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    X_agd_contrast <- X_all[nrow(dat_ipd) + nrow(idat_agd_arm) + 1:nrow(idat_agd_contrast), ]

    # Need to difference .trtb terms - write custom model.matrix?
    abort("Contrast-based AgD not yet supported.")
  } else {
    X_agd_contrast <- NULL
  }

  # Construct RE correlation matrix
  if (trt_effects == "random") {
    dat_all <- dplyr::bind_rows(dat_ipd, dat_agd_arm, dat_agd_contrast)
    if (consistency == "consistency") {
      .RE_cor <- RE_cor(dat_all$.study, dat_all$.trt, type = "reftrt")
      .which_RE <- which_RE(dat_all$.study, dat_all$.trt, type = "reftrt")
    } else if (consistency == "ume") {
      .RE_cor <- RE_cor(dat_all$.study, dat_all$.trt, type = "blshift")
      .which_RE <- which_RE(dat_all$.study, dat_all$.trt, type = "blshift")
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
    n_int = n_int,
    trt_effects = trt_effects,
    RE_cor = .RE_cor,
    which_RE = .which_RE,
    likelihood = likelihood,
    link = link,
    ...,
    prior_intercept = prior_intercept,
    prior_trt = prior_trt,
    prior_het = prior_het,
    prior_reg = prior_reg,
    prior_aux = prior_aux,
    QR = QR,
    adapt_delta = adapt_delta,
    int_thin = int_thin)

  # Create stan_nma object
  out <- list(network = network, stanfit = stanfit)
  class(out) <- "stan_nma"

  return(out)
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
                    n_int,
                    trt_effects = c("fixed", "random"),
                    RE_cor = NULL,
                    which_RE = NULL,
                    likelihood = NULL,
                    link = NULL,
                    ...,
                    prior_intercept = normal(scale = 10),
                    prior_trt = normal(scale = 10),
                    prior_het = half_normal(scale = 5),
                    prior_reg = normal(scale = 10),
                    prior_aux = normal(scale = 5),
                    QR = FALSE,
                    adapt_delta = NULL,
                    int_thin = 100L) {

  if (missing(ipd_x)) ipd_x <- NULL
  if (missing(ipd_y)) ipd_y <- NULL
  if (missing(agd_arm_x)) agd_arm_x <- NULL
  if (missing(agd_arm_y)) agd_arm_y <- NULL
  if (missing(agd_contrast_x)) agd_contrast_x <- NULL
  if (missing(agd_contrast_y)) agd_contrast_y <- NULL

  # Check available x and y
  if (xor(is.null(ipd_x), is.null(ipd_y)))
    abort("`ipd_x` and `ipd_y` should both be present or both NULL.")
  if (xor(is.null(agd_arm_x), is.null(agd_arm_y)))
    abort("`ipd_x` and `ipd_y` should both be present or both NULL.")
  if (xor(is.null(agd_contrast_x), is.null(agd_contrast_y)))
    abort("`ipd_x` and `ipd_y` should both be present or both NULL.")

  has_ipd <- !is.null(ipd_x) && !is.null(ipd_y)
  has_agd_arm <- !is.null(agd_arm_x) && !is.null(agd_arm_y)
  has_agd_contrast <-  !is.null(agd_contrast_x) && !is.null(agd_contrast_y)

  # Check design matrices, outcomes
  if (has_ipd) {
    if (!is.matrix(ipd_x) || !is.numeric(ipd_x))
      abort("`ipd_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(ipd_x, is.numeric)))
      abort("`ipd_y` should be numeric outcome data.")
    if (nrow(ipd_x) != nrow(ipd_y))
      abort("Number of rows in `ipd_x` and `ipd_y` do not match.")
  }
  if (!is.null(agd_arm_x)) {
    if (!is.matrix(agd_arm_x) || !is.numeric(agd_arm_x))
      abort("`agd_arm_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_arm_x, is.numeric)))
      abort("`agd_arm_y` should be numeric outcome data.")
    if (nrow(agd_arm_x) != nrow(agd_arm_y) * n_int)
      abort("Number of rows in `agd_arm_x`, `agd_arm_y`, and `n_int` do not match.")
  }
  if (!is.null(agd_contrast_x)) {
    if (!is.matrix(agd_contrast_x) || !is.numeric(agd_contrast_x))
      abort("`agd_contrast_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_contrast_x, is.numeric)))
      abort("`agd_contrast_y` should be numeric outcome data.")
    if (nrow(agd_contrast_x) != nrow(agd_contrast_y) * n_int)
      abort("Number of rows in `agd_contrast_x`, `agd_contrast_y`, and `n_int` do not match.")
  }

  # Check matching X column names and dimensions if more than one present
  if ((has_ipd &&
       ((has_agd_arm && !identical(colnames(ipd_x), colnames(agd_arm_x)))) ||
       ((has_agd_contrast && !identical(colnames(ipd_x), colnames(agd_contrast_x))))) ||
      (has_agd_arm &&
       (has_agd_contrast && !identical(colnames(agd_arm_x), colnames(agd_contrast_x)))))
    abort("Non-matching columns in *_x matrices.")

  # Check model arguments
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")

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
      trunc(int_thin) != int_thin) abort("`int_thin` should be an integer >= 1.")
  if (!is.numeric(n_int) ||
      length(n_int) > 1 ||
      n_int < 1 ||
      trunc(n_int) != n_int) abort("`n_int` should be an integer >= 1.")

  # Allow n_int = NULL if no AgD
  if (!has_agd_arm && !has_agd_contrast) n_int <- 0

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
  col_trt <- grepl("^\\.trt[^:]+$", x_names)
  col_reg <- !col_study & !col_trt

  n_study <- sum(col_study)
  n_trt <- sum(col_trt) + 1

  get_study <- function(x) which(x == 1)
  get_trt <- function(x, v = 1) if (any(x == v)) which(x == v) + 1 else 1

  if (has_ipd) {
    ipd_s_t <- dplyr::tibble(.study = apply(ipd_x[, col_study], 1, get_study),
                             .trt = apply(ipd_x[, col_trt], 1, get_trt))
    ipd_arm <-  dplyr::group_indices(ipd_s_t, .data$.study, .data$.trt)
    ipd_s_t <- dplyr::distinct(ipd_s_t)
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
    agd_arm_study <- apply(agd_arm_x[aa1, col_study], 1, get_study)
    agd_arm_trt <- apply(agd_arm_x[aa1, col_trt], 1, get_trt)
  } else {
    agd_arm_study <- agd_arm_trt <- numeric()
    ni_agd_arm <- 0
  }

  if (has_agd_contrast) {
    ni_agd_contrast <- nrow(agd_contrast_y)
    ac1 <- 0:(ni_agd_contrast - 1)*n_int + 1
    agd_contrast_study <- apply(agd_contrast_x[ac1, col_study], 1, get_study)
    agd_contrast_trt <- apply(agd_contrast_x[ac1, col_trt], 1, get_trt)
    agd_contrast_trt_b <- apply(agd_contrast_x[ac1, col_trt], 1, get_trt, v = -1)
  } else {
    agd_contrast_study <- agd_contrast_trt <- agd_contrast_trt_b <- numeric()
    ni_agd_contrast <- 0
  }

  # Set up random effects
  if (trt_effects == "random") {
    narm <- narm_ipd + ni_agd_arm + ni_agd_contrast
    if (!is.null(which_RE)) {
      if (!is.numeric(which_RE) ||
          trunc(which_RE) != which_RE ||
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
    # ns_agd_contrast = length(unique(agd_contrast_study)),
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
    agd_contrast_trt = agd_contrast_trt,
    agd_contrast_trt_b = agd_contrast_trt_b,
    agd_contrast_y = if (has_agd_contrast) agd_contrast_y$.y else numeric(),
    agd_contrast_se = if (has_agd_contrast) agd_contrast_y$.se else numeric(),
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
    R_inv = if (QR) X_all_R_inv else matrix(0, 0, 0)
    )

  # Add priors
  standat <- purrr::list_modify(standat,
    !!! prior_standat(prior_intercept, "prior_intercept",
                      valid = c("Normal", "Cauchy", "Student t")),
    !!! prior_standat(prior_trt, "prior_trt",
                      valid = c("Normal", "Cauchy", "Student t")),
    !!! prior_standat(prior_reg, "prior_reg",
                      valid = c("Normal", "Cauchy", "Student t")),
    !!! prior_standat(prior_het, "prior_het",
                      valid = c("Normal", "half-Normal",
                                "Cauchy",  "half-Cauchy",
                                "Student t", "half-Student t"))
    )

  # Standard pars to monitor
  pars <- c("mu", "beta", "gamma",
            "log_lik", "resdev", "lp__")

  # Monitor heterogeneity SD and study deltas if RE model
  if (trt_effects == "random") {
    pars <- c(pars, "tau", "delta")
  }
  # Monitor cumulative integration error if using numerical integration
  if (n_int > 1) {
    pars <- c(pars, "theta_bar_cum")
  }

  # Call Stan model for given likelihood

  # -- Normal likelihood
  if (likelihood == "normal") {

    standat <- purrr::list_modify(standat,
    # Add outcomes
      ipd_y = ipd_y$.y,
      agd_arm_y = agd_arm_y$.y, agd_arm_se = agd_arm_y$.se,

    # Add prior for auxilliary parameter - individual-level variance
      !!! prior_standat(prior_het, "prior_aux",
                        valid = c("Normal", "half-Normal",
                                  "Cauchy",  "half-Cauchy",
                                  "Student t", "half-Student t")),

    # Specify link
    link = switch(link, identity = 1, log = 2)
    )

    stanfit <- rstan::sampling(stanmodels$normal, data = standat,
                               pars = c(pars, "sigma"), ...)


  # -- Bernoulli/binomial likelihood (one parameter)
  } else if (likelihood %in% c("bernoulli", "binomial")) {

    standat <- purrr::list_modify(standat,
    # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2)
    )

    stanfit <- rstan::sampling(stanmodels$binomial_1par, data = standat,
                               pars = pars, ...)

  # -- Bernoulli/binomial likelihood (two parameter)
  } else if (likelihood %in% c("bernoulli2", "binomial2")) {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2)
    )


    stanfit <- rstan::sampling(stanmodels$binomial_2par, data = standat,
                               pars = c(pars, "theta2_bar_cum"), ...)

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

    stanfit <- rstan::sampling(stanmodels$poisson, data = standat,
                               pars = pars, ...)

  } else {
    abort(glue::glue('"{likelihood}" likelihood not supported.'))
  }

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
#' @param type Character string, whether to generate RE structure under the
#'   "reference treatment" paramterisation, or the "baseline shift"
#'   parameterisation.
#'
#' @return
#' @export
#' @aliases RE_cor
#' @rdname random_effects
#'
#' @examples
#' RE_cor(smoking$studyn, smoking$trtn)
#' RE_cor(smoking$studyn, smoking$trtn, type = "blshift")
RE_cor <- function(study, trt, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (length(study) != length(trt)) abort("`study` and `trt` must be the same length.")
  type <- rlang::arg_match(type)


  if (type == "reftrt") {
    reftrt <- levels(trt)[1]
    nonref <- trt != reftrt
    nRE <- sum(nonref)  # RE for each non ref trt arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[nonref]
    trt <- trt[nonref]
  } else if (type == "blshift") {
    bl <- !duplicated(study)
    nRE <- sum(!bl)  # RE for each non baseline arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[!bl]
    trt <- trt[!bl]
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
#' which_RE(smoking$studyn, smoking$trtn)
#' which_RE(smoking$studyn, smoking$trtn, type = "blshift")
which_RE <- function(study, trt, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (length(study) != length(trt)) abort("`study` and `trt` must be the same length.")
  type <- rlang::arg_match(type)

  n_i <- length(study)
  id <- rep(0, n_i)

  if (type == "reftrt") {
    reftrt <- levels(trt)[1]
    trt_nonref <- trt != reftrt
    id[trt_nonref] <- 1:sum(trt_nonref)
  } else if (type == "blshift") {
    non_bl <- duplicated(study)
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
                      continuous = "normal")

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
                     glue::glue_collapse(dQuote(valid_lhood, FALSE), sep = ", ", last = " or "),
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
                     bernoulli = c("logit", "probit"),
                     bernoulli2 = c("logit", "probit"),
                     binomial = c("logit", "probit"),
                     binomial2 = c("logit", "probit"),
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
                glue::collapse(valid, sep = ", ", last = ", or ")))

  type <- switch(dist,
                 Normal = , `half-Normal` = 1,
                 Cauchy = , `half-Cauchy` = 2,
                 `Student t` = , `half-Student t` = 3,
                 Exponential = 4)

  out <- purrr::list_modify(c(x), type = type, dist = purrr::zap())
  # Set unnecessary (NA) parameters to zero. These will be ignored by Stan, but
  # need to pass rstan checks
  out[is.na(out)] <- 0
  names(out) <- paste0(par, "_", names(out))
  return(out)
}
