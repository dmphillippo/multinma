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
    # rlang::with_handlers(
    out$agd_arm <-
      dplyr::bind_cols(
        network$agd_arm,
        purrr::pmap_dfc(list(x_int_names, ds, u_cor_l),
                        ~ rowwise(network$agd_arm) %>%
                          transmute(!! ..1 := list(rlang::eval_tidy(rlang::call2(..2$qfun, p = ..3, !!! ..2$args)))))
      )#,
    # error = abort,
    # warning = rlang::calling(~{warn(.); rlang::cnd_muffle(.)})
    # )

    # Check valid values produced
    invalid_rows <- out$agd_arm %>%
      dplyr::mutate_at(x_int_names,
                       .funs = ~purrr::map_lgl(.,
                                 ~any(is.na(.) || is.infinite(.) ||
                                      is.null(.) || is.nan(.)))
                       ) %>%
      dplyr::filter_at(x_int_names, dplyr::any_vars(.))

    if (nrow(invalid_rows) > 0) {
      abort(
        glue::glue("Invalid integration points were generated (either NA, NaN, Inf, or NULL).\n",
                   "Check the input parameters for the following (study, treatment):\n",
                   glue::glue_collapse(glue::glue(" {invalid_rows$.study}, {invalid_rows$.trt}"),
                                       sep = "\n"))
      )
    }
  }

  if (has_agd_contrast(network)) {
    # rlang::with_handlers(
    out$agd_contrast <-
      dplyr::bind_cols(
        network$agd_contrast,
        purrr::pmap_dfc(list(x_int_names, ds, u_cor_l),
                        ~ rowwise(network$agd_contrast) %>%
                          transmute(!! ..1 := list(rlang::eval_tidy(rlang::call2(..2$qfun, p = ..3, !!! ..2$args)))))
      )#,
    # error = abort,
    # warning = rlang::calling(~{warn(.); rlang::cnd_muffle(.)})
    # )

    # Check valid values produced
    invalid_rows <- out$agd_contrast %>%
      dplyr::mutate_at(x_int_names,
                       .funs = ~purrr::map_lgl(.,
                                               ~any(is.na(.) || is.infinite(.) ||
                                                      is.null(.) || is.nan(.)))
      ) %>%
      dplyr::filter_at(x_int_names, dplyr::any_vars(.))

    if (nrow(invalid_rows) > 0) {
      abort(
        glue::glue("Invalid integration points were generated (either NA, NaN, Inf, or NULL).\n",
                   "Check the input parameters for the following (study, comparison):\n",
                   glue::glue_collapse(glue::glue(" {invalid_rows$.study}, {invalid_rows$.trt} vs. {invalid_rows$.trt_b}"),
                                       sep = "\n"))
      )
    }
  }

  # Set as mlnmr_data class
  out$n_int <- n_int
  out$int_names <- x_names
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

#' The Bernoulli Distribution
#'
#' The quantile function `qbern` for a Bernoulli distribution, with success
#' probability `prob`. This is equivalent to `qbinom(p, 1, prob)`.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param prob probability of success
#' @param lower.tail,log.p,log see [stats::Binomial]
#'
#' @export
#' @rdname Bernoulli
#' @aliases qbern
qbern <- function(p, prob, lower.tail = TRUE, log.p = FALSE) {
  return(qbinom(p, size = 1, prob = prob, lower.tail = lower.tail, log.p = log.p))
}

#' @export
#' @rdname Bernoulli
#' @aliases pbern
pbern <- function(q, prob, lower.tail = TRUE, log.p = FALSE) {
  return(pbinom(q, size = 1, prob = prob, lower.tail = lower.tail, log.p = log.p))
}

#' @export
#' @rdname Bernoulli
#' @aliases dbern
dbern <- function(x, prob, log = FALSE) {
  return(dbinom(x, size = 1, prob = prob, log = log))
}

#' The Gamma distribution
#'
#' We provide convenient extensions of the `[dpq]gamma` functions, which allow
#' the distribution to be specified in terms of its mean and standard deviation,
#' instead of shape and rate/scale.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param shape,rate,scale,log,lower.tail,log.p see [stats::GammaDist]
#' @param mean,sd mean and standard deviation, overriding `shape` and
#'   `rate` or `scale` if specified
#'
#' @return
#' @export
#' @rdname GammaDist
#' @aliases qgamma
qgamma <- function(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE,
                   log.p = FALSE, mean, sd) {
  if (!missing(mean) && !missing(sd))
    return(stats::qgamma(p, shape = (mean / sd)^2, rate = mean / sd^2,
                        lower.tail = lower.tail, log.p = log.p))
  else
    return(stats::qgamma(p, shape = shape, rate = rate,
                        lower.tail = lower.tail, log.p = log.p))
}

#' @export
#' @rdname GammaDist
#' @aliases dgamma
dgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE, mean, sd) {
  if (!missing(mean) && !missing(sd))
    return(stats::dgamma(x, shape = (mean / sd)^2, rate = mean / sd^2, log = log))
  else
    return(stats::dgamma(x, shape = shape, rate = rate, log = log))
}

#' @export
#' @rdname GammaDist
#' @aliases pgamma
pgamma <- function(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE,
                   log.p = FALSE, mean, sd) {
  if (!missing(mean) && !missing(sd))
    return(stats::pgamma(q, shape = (mean / sd)^2, rate = mean / sd^2,
                         lower.tail = lower.tail, log.p = log.p))
  else
    return(stats::pgamma(q, shape = shape, rate = rate,
                         lower.tail = lower.tail, log.p = log.p))
}


#' The logit Normal distribution
#'
#' We provide convenient extensions of the `[dpq]logitnorm` functions in the
#' package [logitnorm], which allow the distribution to be specified in terms of
#' its mean and standard deviation, instead of its logit-mean and logit-sd.
#'
#' @param p,x vector of quantiles
#' @param q vector of probabilities
#' @param mu,sigma see [logitnorm]
#' @param ...
#' @param mean,sd mean and standard deviation, overriding `mu` and `sigma` if
#'   specified
#'
#' @return
#' @export
#' @rdname logitNormal
#' @aliases qlogitnorm
qlogitnorm <- function(p, mu = 0, sigma = 1, ..., mean, sd){
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::qlogitnorm(p, pars[["mu"]], pars[["sigma"]]))
}

#' @export
#' @rdname logitNormal
#' @aliases dlogitnorm
dlogitnorm <- function(x, mu = 0, sigma = 1, ..., mean, sd) {
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::dlogitnorm(x, pars[["mu"]], pars[["sigma"]]))
}

#' @export
#' @rdname logitNormal
#' @aliases plogitnorm
plogitnorm <- function(q, mu = 0, sigma = 1, ..., mean, sd) {
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::plogitnorm(q, pars[["mu"]], pars[["sigma"]]))
}

# Internal functions for *logitnorm()
.lndiff <- function(est, m, s){
  x <- logitnorm::momentsLogitnorm(est[1], est[2])
  return((x[1] - m)^2 + (sqrt(x[2]) - s)^2)
}

.lnopt <- function(m, s) {
  opt <- optim(c(m, s), .lndiff, m = m, s = s)

  if (opt$convergence != 0) {
    warn("Optimisation did not converge, NAs produced.")
    return(c(mu = NA, sigma = NA))
  } else {
    return(c(mu = opt$par[1], sigma = opt$par[2]))
  }
}

# Estimate mu and sigma parameters of logit-normal from sample mean and sd
pars_logitnorm <- function(m, s) {
  if (length(m) != length(s)) abort("Parameters not same length.")
  if (any(m > 1 | m < 0)) abort("Mean not in [0,1]. Have you rescaled?")

  return(as.data.frame(do.call(rbind, mapply(.lnopt, m, s, SIMPLIFY = FALSE))))

}
