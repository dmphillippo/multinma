#' Add numerical integration points to aggregate data
#'
#' The `add_integration()` generic creates numerical integration points using a
#' Gaussian copula approach, as described in
#' \insertCite{methods_paper;textual}{multinma}. Methods are available for
#' networks stored in `nma_data` objects, and for data frames. The function
#' `unnest_integration()` unnests integration points stored in a data frame, to
#' aid plotting or other exploration.
#'
#' @param x An `nma_data` object, as created by the `set_*()` functions or
#'   `combine_network()`, or data frame
#' @param ... Distributions for covariates, see "Details"
#' @param cor Correlation matrix to use for generating the integration points.
#'   By default, this takes a weighted correlation matrix from all IPD studies.
#'   Rows and columns should match the order of covariates specified in `...`.
#' @param n_int Number of integration points to generate, default 1000
#' @param int_args A named list of arguments to pass to
#'   \code{\link[randtoolbox:quasiRNG]{sobol()}}
#'
#' @return For the `nma_data` method, an object of class [nma_data]. For the
#'   `data.frame` method, the input data frame is returned (as a [tibble]) with
#'   an added column for each covariate (prefixed with ".int_"), containing the
#'   numerical integration points nested as length-`n_int` vectors within each
#'   row. For `unnest_integration()`, a data frame with integration points
#'   unnested.
#' @export
#'
#' @details The arguments passed to `...` specify distributions for the
#'   covariates. Argument names specify the name of the covariate, which should
#'   match a covariate name in the IPD (if IPD are present). The required
#'   marginal distribution is then specified using the function [distr()].
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples ## Plaque psoriasis ML-NMR - network setup and adding integration points
#' @template ex_plaque_psoriasis_network
#' @template ex_plaque_psoriasis_integration
#' @examples
#'
#' ## Adding integration points to a data frame, e.g. for prediction
#' # Define a data frame of covariate summaries
#' new_agd_int <- data.frame(
#'   bsa_mean = 0.6,
#'   bsa_sd = 0.3,
#'   prevsys = 0.1,
#'   psa = 0.2,
#'   weight_mean = 10,
#'   weight_sd = 1,
#'   durnpso_mean = 3,
#'   durnpso_sd = 1)
#'
#' # Adding integration points, using the weighted average correlation matrix
#' # computed for the plaque psoriasis network
#' new_agd_int <- add_integration(new_agd_int,
#'   durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
#'   prevsys = distr(qbern, prob = prevsys),
#'   bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
#'   weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
#'   psa = distr(qbern, prob = psa),
#'   cor = pso_net$int_cor,
#'   n_int = 1000)
#'
#' new_agd_int
#'
add_integration <- function(x, ...) {
  UseMethod("add_integration")
}

#' @export
#' @rdname add_integration
add_integration.default <- function(x, ...) {
  abort(glue::glue("No add_integration method defined for class '{class(x)}'."))
}

#' @export
#' @rdname add_integration
add_integration.data.frame <- function(x, ...,
                                       cor = NULL, n_int = 1000L, int_args = list()) {

  x <- tibble::as_tibble(x)

  # Check n_int
  if (!rlang::is_scalar_integerish(n_int) || n_int <= 0) {
    abort("`n_int` should be a positive integer.")
  }

  # Check int_args
  if (!is.list(int_args) || (!rlang::is_empty(int_args) && !rlang::is_named(int_args))) {
    abort("`int_args` should be a named list.")
  }

  # Check covariate arguments
  ds <- list(...)

  if (length(ds) == 0) {
    abort("No covariate distributions specified. Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  if (any(purrr::map_lgl(ds, ~!inherits(., "distr"))) || !rlang::is_named(ds)) {
    abort("Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  x_names <- names(ds)  # Covariate names
  nx <- length(ds)      # Number of covariates

  # Check cor
  if (nx > 1) {
    if (!is.null(cor)) {
      tryCatch(cor <- as.matrix(cor),
               error = function(e) abort("`cor` should be a correlation matrix or NULL"))

      if (!is.numeric(cor) ||
          !isSymmetric(cor) ||
          !isTRUE(all.equal(diag(cor), rep(1, nrow(cor)), check.names = FALSE)) ||
          !all(eigen(cor, symmetric = TRUE)$values > 0))
        abort("`cor` should be a correlation matrix or NULL")

      if (ncol(cor) != nx)
        abort("Dimensions of correlation matrix `cor` and number of covariates specified in `...` do not match.")
    } else {
      abort("Specify a correlation matrix using the `cor` argument.")
    }
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
  out <- x
  x_int_names <- paste0(".int_", x_names)

  # If integration points already present, remove and warn
  int_cols <- stringr::str_detect(colnames(x), "^\\.int_")
  if (any(int_cols))
    warn("Replacing integration points already present in data frame.", "int_col_present")

  # Bind in generated integration points
  out <-
    dplyr::bind_cols(
      x[, !int_cols],
      purrr::pmap_dfc(list(x_int_names, ds, u_cor_l),
                      ~ dplyr::rowwise(x) %>%
                        dplyr::transmute(!! ..1 := list(rlang::eval_tidy(rlang::call2(..2$qfun, p = ..3, !!! ..2$args)))))
    )

  # Check valid values produced
  invalid_rows <- out %>%
    dplyr::mutate_at(x_int_names,
                     .funs = ~purrr::map_lgl(.,
                                             ~any(is.na(.) | is.infinite(.) |
                                                    is.null(.) | is.nan(.)))
    ) %>%
    dplyr::filter_at(x_int_names, dplyr::any_vars(.))

  if (nrow(invalid_rows) > 0)
    abort("Invalid integration points were generated (either NA, NaN, Inf, or NULL).",
          "invalid_int_generated",
          invalid_rows = invalid_rows)

  class(out) <- c("integration_tbl", class(out))
  return(out)
}

#' @export
#' @rdname add_integration
add_integration.nma_data <- function(x, ...,
                                     cor = NULL, n_int = 1000L, int_args = list()) {

  network <- x

  # Checks for n_int and int_args performed in data.frame method

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

  # Covariate arguments
  ds <- list(...)

  if (length(ds) == 0) {
    abort("No covariate distributions specified. Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  if (any(purrr::map_lgl(ds, ~!inherits(., "distr"))) || !rlang::is_named(ds)) {
    abort("Covariate distributions should be specified as named arguments using the function `distr`.")
  }

  x_names <- names(ds)
  nx <- length(ds)

  # Check cor
  if (!is.null(cor)) {
    tryCatch(cor <- as.matrix(cor),
             error = function(e) abort("`cor` should be a correlation matrix or NULL"))

    if (!is.numeric(cor) ||
        !isSymmetric(cor) ||
        !isTRUE(all.equal(diag(cor), rep(1, nrow(cor)), check.names = FALSE)) ||
        !all(eigen(cor, symmetric = TRUE)$values > 0)) {
      abort("`cor` should be a correlation matrix or NULL")
    }
  } else if (nx > 1) {
    if (!has_ipd(network)) abort("Specify a correlation matrix using the `cor` argument, or provide IPD studies in the network.")
  }

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

  out <- network

  # Establish handlers for more specific errors / warnings when possible
  int_col_present <- function(wrn) {
    warn("Replacing integration points already present in network.")
    rlang::cnd_muffle(wrn)
  }

  invalid_int_generated <- function(err) {
    invalid_rows <- err$invalid_rows
    abort(
      glue::glue("Invalid integration points were generated (either NA, NaN, Inf, or NULL).\n",
                 "Check the input parameters for the following (study, treatment):\n",
                 glue::glue_collapse(glue::glue(" {invalid_rows$.study}, {invalid_rows$.trt}"),
                                     sep = "\n"))
    )
  }

  if (has_agd_arm(network)) {
    out$agd_arm <- rlang::with_handlers(
      add_integration.data.frame(network$agd_arm, ...,
                                 cor = cor, n_int = n_int, int_args = int_args),
      int_col_present = rlang::calling(int_col_present),
      invalid_int_generated = invalid_int_generated)
  }

  if (has_agd_contrast(network)) {
    out$agd_contrast <- rlang::with_handlers(
      add_integration.data.frame(network$agd_contrast, ...,
                                 cor = cor, n_int = n_int, int_args = int_args),
      int_col_present = rlang::calling(int_col_present),
      invalid_int_generated = invalid_int_generated)
  }

  # Set as mlnmr_data class
  out$n_int <- n_int
  out$int_names <- x_names

  if (nx == 1) cor <- matrix(1)
  colnames(cor) <- x_names
  rownames(cor) <- x_names
  out$int_cor <- cor

  class(out) <- c("mlnmr_data", "nma_data")
  return(out)
}


#' @param data Data frame with nested integration points, stored in list
#'   columns as `.int_<variable name>`
#'
#' @export
#' @rdname add_integration
unnest_integration <- function(data) {
  if (!inherits(data, "data.frame")) abort("Expecting a data frame.")

  # Get integration variable names
  data_names <- colnames(data)
  x_int_names <- stringr::str_subset(data_names, "^\\.int_")
  x_names <- stringr::str_remove(x_int_names, "^\\.int_")

  if (length(x_int_names) == 0)
    abort("No integration points present in data frame.")

  name_conflicts <- rlang::has_name(data, x_names)
  if (any(name_conflicts)) {
    warn(glue::glue("Columns will be overwritten when unnesting integration points: ",
                    glue::glue_collapse(x_names[name_conflicts], sep = ", ")),
         "unnest_name_conflict")
  }

  out <- data %>%
    dplyr::select(-dplyr::one_of(x_names[name_conflicts])) %>%
    dplyr::rename_at(x_int_names, ~stringr::str_remove(., "^\\.int_")) %>%
    {if (getNamespaceVersion("tidyr") < "1.0.0") {
      tidyr::unnest(., !!! rlang::syms(x_names))
    } else {
      tidyr::unnest(., x_names)
    }}

  return(out)
}

#' Internal version of unnest_integration, doesn't throw warnings when
#' overwriting columns
#' @noRd
.unnest_integration <- function(data) {
  ignore <- function(wrn) rlang::cnd_muffle(wrn)
  out <- rlang::with_handlers(unnest_integration(data),
                              unnest_name_conflict = rlang::calling(ignore))
  return(out)
}

#' Specify a general marginal distribution
#'
#' `distr()` is used within the function [add_integration()] to specify marginal
#' distributions for the covariates, via a corresponding inverse CDF. It is also
#' used in [predict.stan_nma()] to specify a distribution for the baseline
#' response (intercept) when predicting absolute outcomes.
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
#' @seealso [add_integration()] where `distr()` is used to specify marginal
#'   distributions for covariates to integrate over, and [predict.stan_nma()]
#'   where `distr()` is used to specify a distribution on the baseline response.
#' @examples
#' ## Specifying marginal distributions for integration
#'
#' df <- data.frame(x1_mean = 2, x1_sd = 0.5, x2 = 0.8)
#'
#' # Distribution parameters are evaluated in the context of the data frame
#' add_integration(df,
#'                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
#'                 x2 = distr(qbern, prob = x2),
#'                 cor = diag(2))
#'
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
#' package \link[logitnorm:logitnorm-package]{logitnorm}, which allow the
#' distribution to be specified in terms of its mean and standard deviation,
#' instead of its logit-mean and logit-sd.
#'
#' @param p,x vector of quantiles
#' @param q vector of probabilities
#' @param mu,sigma,... see \code{\link[logitnorm:logitnorm-package]{logitnorm}}
#' @param mean,sd mean and standard deviation, overriding `mu` and `sigma` if
#'   specified
#'
#' @export
#' @rdname logitNormal
#' @aliases qlogitnorm
qlogitnorm <- function(p, mu = 0, sigma = 1, ..., mean, sd){
  require_pkg("logitnorm")
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::qlogitnorm(p, pars[["mu"]], pars[["sigma"]], ...))
}

#' @export
#' @rdname logitNormal
#' @aliases dlogitnorm
dlogitnorm <- function(x, mu = 0, sigma = 1, ..., mean, sd) {
  require_pkg("logitnorm")
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::dlogitnorm(x, pars[["mu"]], pars[["sigma"]], ...))
}

#' @export
#' @rdname logitNormal
#' @aliases plogitnorm
plogitnorm <- function(q, mu = 0, sigma = 1, ..., mean, sd) {
  require_pkg("logitnorm")
  if (!missing(mean) && !missing(sd)) pars <- pars_logitnorm(mean, sd)
  else pars <- list(mu = mu, sigma = sigma)
  return(logitnorm::plogitnorm(q, pars[["mu"]], pars[["sigma"]], ...))
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
