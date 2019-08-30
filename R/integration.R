#' Add numerical integration points to aggregate data
#'
#' This function creates numerical integration points using a Gaussian copula
#' approach, as described in \insertCite{methods_paper}{multinma}.
#'
#' @param network an `nma_data` object, as created by the `set_*()` functions or
#'   `combine_network()`
#' @param ... distributions for covariates, see "Details"
#' @param cor correlation matrix to use for generating the integration points.
#'   By default, this takes a weighted correlation matrix from all IPD studies.
#' @param n_int number of integration points to generate, default 100
#' @param int_args a named list of arguments to pass to [randtoolbox::sobol()]
#'
#' @return An object of class [nma_data].
#' @export
#'
#' @details The arguments passed to `...` specify distributions for the
#'   covariates. Argument names specify the name of the covariate, which should
#'   match a covariate name in the IPD (if IPD are present). The required
#'   marginal distribution is then specified using the function [distr()].
#'
#' @examples
add_integration <- function(network, ..., cor = NULL, n_int = 100L, int_args = list()) {
  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the `set_*` or `combine_network` functions.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  if (!has_agd_arm(network) && !has_agd_contrast(network)) {
    abort("No aggregate data found in network.")
  }

  # Check n_int
  if (!is.numeric(n_int) || n_int != trunc(n_int) || n_int <= 0 || length(n_int) > 1) {
    abort("`n_int` should be a positive integer.")
  }

  # Check int_args
  if (!is.list(int_args) || (!rlang::is_empty(int_args) && !rlang::is_named(int_args))) {
    abort("`int_args` should be a named list.")
  }

  # Check cor
  if (!is.null(cor)) {
    tryCatch(cor <- as.matrix(cor),
             error = function(e) abort("cor should be a correlation matrix or NULL"))

    if (!is.numeric(cor) ||
        !isSymmetric(cor) ||
        !isTRUE(all.equal(diag(cor), rep(1, nrow(cor)))) ||
        !all(eigen(cor, symmetric = TRUE)$values > 0)) {
      abort("cor should be a correlation matrix or NULL")
    }
  } else {
    if (!has_ipd(network)) abort("Specify a correlation matrix using the `cor` argument, or provide IPD studies in the network.")
  }

  # Check covariate arguments
  ds <- list(...)

  if (length(ds) == 0) {
    abort("No covariate distributions specified. Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  if (any(purrr::map_lgl(ds, ~!inherits(., "distr"))) || !rlang::is_named(ds)) {
    abort("Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  x_names <- names(ds)
  nx <- length(ds)

  # If IPD is provided, check that covariate names match
  if (has_ipd(network) && any(! x_names %in% colnames(network$ipd))) {
    abort(paste0("Covariate name(s) not found in IPD: ",
                 paste0(setdiff(x_names, colnames(network$ipd)), collapse = ", ")))
  }

  # Use weighted average correlation matrix from IPD, if cor = NULL
  if (is.null(cor) && nx > 1) {
    inform("Using weighted average correlation matrix computed from IPD studies.")


    # Check for any missing covariates
    if (!all(complete.cases(dplyr::select(network$ipd, !! x_names)))) {
      warn("Missing values found for some covariates in IPD. Calculating correlations using complete cases.")
    }

    # Weighted z-transform average
    ipd_cors <-
      network$ipd %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::group_modify(~tibble::tibble(
        w = nrow(.) - 3,
        r = list(cor(dplyr::select(., !! x_names), method = "spearman", use = "complete.obs"))
        )) %>%
      dplyr::mutate(z = purrr::map2(.data$w, .data$r, ~.x * log((1 + .y) / (1 - .y)) / 2))

    w_z <- Reduce(`+`, ipd_cors$z) / sum(ipd_cors$w)
    ipd_cor <- (exp(2 * w_z) - 1) / (exp(2 * w_z) + 1)

    diag(ipd_cor) <- 1

    cor <- ipd_cor
  }

  # Generate Sobol points
  u <- do.call(randtoolbox::sobol, purrr::list_modify(list(n = n_int, dim = nx), !!! int_args))

  # Correlate Sobol points with Gaussian copula
  if (nx > 1) {
    cop <- copula::normalCopula(copula::P2p(cor), dim = nx, dispstr = "un")
    u_cor <- copula::cCopula(u, copula = cop, inverse = TRUE)
    # columns to list
    u_cor_l <- purrr::array_branch(u_cor, 2)
  } else {
    u_cor <- u
    u_cor_l <- list(u_cor)
  }

  # Use inverse CDFs to get integration points for specified marginals
  out <- network
  x_int_names <- paste0(".int_", x_names)

  if (has_agd_arm(network)) {
    out$agd_arm <-
      dplyr::bind_cols(
        network$agd_arm,
        purrr::pmap_dfc(list(x_int_names, ds, u_cor_l),
                        ~ rowwise(network$agd_arm) %>%
                          transmute(!! ..1 := list(rlang::eval_tidy(rlang::call2(..2$qfun, p = ..3, !!! ..2$args)))))
      )
  }

  if (has_agd_contrast(network)) {
    out$agd_contrast <-
      dplyr::bind_cols(
        network$agd_contrast,
        purrr::pmap_dfc(list(x_int_names, ds, u_cor_l),
                        ~ rowwise(network$agd_contrast) %>%
                          transmute(!! ..1 := list(rlang::eval_tidy(rlang::call2(..2$qfun, p = ..3, !!! ..2$args)))))
      )
  }

  # Set as mlnmr_data class
  out$n_int <- n_int
  class(out) <- c("mlnmr_data", "nma_data")
  return(out)
}


#' Specify a general marginal distribution
#'
#' `distr()` is used within the function [add_integration()] to specify marginal
#' distributions for the covariates, via a corresponding inverse CDF.
#'
#' @param qfun an inverse CDF, either as a function name or a string
#' @param ... parameters of the distribution as arguments to `qfun`, these will
#'   be quoted and evaluated later in the context of the aggregate data sources
#'
#' @return An object of class [distr].
#' @export
#'
#' @details The function `qfun` should have a formal argument called `p`. This
#'   restriction serves as a crude check for inverse CDFs (e.g. an error will be
#'   given if `dnorm` is used instead of `qnorm`). If a user-written CDF is
#'   supplied, it must have an argument `p` which takes a vector of
#'   probabilities.
#'
#' @examples
distr <- function(qfun, ...) {
  d <- list(qfun = match.fun(qfun),
            args = rlang::enquos(...))
  if (! "p" %in% rlang::fn_fmls_names(qfun)) {
    abort("`qfun` should be an inverse CDF function (for example `qnorm`, `qgamma`, `qbinom`, ...) but does not appear to be (no formal argument `p`)")
  }
  class(d) <- "distr"
  return(d)
}
