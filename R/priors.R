#' Prior distributions
#'
#' These functions are used to specify prior distributions for the model parameters.
#'
#' @param location Prior location. Typically prior mean (see details).
#' @param scale Prior scale. Typically prior standard deviation (see details).
#' @param df Prior degrees of freedom.
#' @param rate Prior rate.
# #' @param lower,upper Lower and upper bounds for a uniform prior distribution.
#'
#' @rdname priors
#' @name priors
#' @aliases priors normal
#'
#' @details The `location` and `scale` parameters are typically the prior mean
#'   and standard deviation, with the following exceptions:
#'   \itemize{
#'   \item For the Cauchy distribution `location` is the prior median and
#'   `scale` is the prior scale.
#'   \item For the log-Normal distribution, `location` and `scale` are the prior
#'   mean and standard deviation of the logarithm.
#'   }
#'
#' ## Compatibility with model parameters
#' The following table summarises which prior distributions may be used with
#' which model parameters. Essentially, priors that take only non-negative
#' values (e.g. half-Normal) may only be used for non-negative parameters
#' (heterogeneity SD/variance/precision, and any auxiliary parameter). If a
#' real-valued prior distribution is specified for a non-negative parameter, it
#' will be truncated at 0 to be non-negative.
#'
#' |       | \strong{Intercept} `prior_intercept` | \strong{Treatment effects} `prior_trt` | \strong{Heterogeneity} `prior_het` | \strong{Regression coefficients} `prior_reg` | \strong{Auxiliary parameter} `prior_aux` |
#' | ----- | :---: | :---: | :---: | :---: | :---: |
#' | \strong{Normal} `normal()` | Yes | Yes | Yes | Yes | Yes |
#' | \strong{half-Normal} `half_normal()` | - | - | Yes | - | Yes |
#' | \strong{log-Normal} `log_normal()` | - | - | Yes | - | Yes |
#' | \strong{Cauchy }`cauchy()` | Yes | Yes | Yes | Yes | Yes |
#' | \strong{half-Cauchy} `half_cauchy()` | - | - | Yes | - | Yes |
#' | \strong{Student t} `student_t()` | Yes | Yes | Yes | Yes | Yes |
#' | \strong{half-Student t} `half_student_t()` | - | - | Yes | - | Yes |
#' | \strong{log-Student t} `log_student_t()` | - | - | Yes | - | Yes |
#' | \strong{Exponential} `exponential()` | - | - | Yes | - | Yes |
#' | \strong{Flat} `flat()` | Yes | Yes | Yes | Yes | Yes |
#'
#' The `flat()` prior is a special case where no prior information is added to
#' the model, resulting in an implicit flat uniform prior distribution over the
#' entire support for a parameter. This will be an improper prior if the
#' parameter is unbounded, and is not generally advised. See the
#' [Stan user's guide](https://mc-stan.org/docs/stan-users-guide/some-differences-in-the-statistical-models-that-are-allowed.html)
#' for more details.
#'
#' @return Object of class [nma_prior].
#' @export
#'
#' @seealso [summary.nma_prior()] for summarising details of prior
#'   distributions. [plot_prior_posterior()] for plots comparing the prior and
#'   posterior distributions of model parameters.
normal <- function(location = 0, scale) {
  check_prior_location(location)
  check_prior_scale(scale)

  return(new_nma_prior("Normal", location = location, scale = scale))
}

#' @rdname priors
#' @export
half_normal <- function(scale) {
  check_prior_scale(scale)

  return(new_nma_prior("half-Normal", location = 0, scale = scale))
}

#' @rdname priors
#' @export
log_normal <- function(location, scale) {
  check_prior_location(location)
  check_prior_scale(scale)

  return(new_nma_prior("log-Normal", location = location, scale = scale))
}

#' @rdname priors
#' @export
cauchy <- function(location = 0, scale) {
  check_prior_location(location, type = "location (median)")
  check_prior_scale(scale, type = "scale")

  return(new_nma_prior("Cauchy", location = location, scale = scale))
}

#' @rdname priors
#' @export
half_cauchy <- function(scale) {
  check_prior_scale(scale, type = "scale")

  return(new_nma_prior("half-Cauchy", location = 0, scale = scale))
}

#' @rdname priors
#' @export
student_t <- function(location = 0, scale, df) {
  check_prior_location(location)
  check_prior_scale(scale)
  check_prior_scale(df, type = "degrees of freedom")

  return(new_nma_prior("Student t", location = location, scale = scale, df = df))
}

#' @rdname priors
#' @export
half_student_t <- function(scale, df) {
  check_prior_scale(scale)
  check_prior_scale(df, type = "degrees of freedom")

  return(new_nma_prior("half-Student t", location = 0, scale = scale, df = df))
}

#' @rdname priors
#' @export
log_student_t <- function(location, scale, df) {
  check_prior_location(location)
  check_prior_scale(scale)
  check_prior_scale(df, type = "degrees of freedom")

  return(new_nma_prior("log-Student t", location = location, scale = scale, df = df))
}

#' @rdname priors
#' @export
exponential <- function(scale = 1/rate, rate = 1/scale) {
  if (missing(scale) && missing(rate))
    abort("Missing argument. Specify either `rate` or `scale`.")
  if (!missing(rate) && !missing(scale))
    warn("Both `rate` and `scale` provided, only `scale` will be used")

  check_prior_scale(scale, type = "scale or rate")

  return(new_nma_prior("Exponential", scale = scale))
}

#' @rdname priors
#' @export
flat <- function() {
  return(new_nma_prior("flat (implicit)"))
}

# #' @rdname priors
# #' @export
# uniform <- function(lower, upper) {
#   check_prior_location(lower, type = "lower limit")
#   check_prior_location(lower, type = "upper limit")
#
#   return(new_nma_prior("Uniform", lower = lower, upper = upper))
# }

# Check functions
check_prior_location <- function(x, type = "location (mean)") {
  if (!(.is_default(x) && rlang::is_empty(x))) {
    if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
    if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
  }
}

check_prior_scale <- function(x, type = "scale (standard deviation)") {
  if (!(.is_default(x) && rlang::is_empty(x))) {
    if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
    if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
    if (x <= 0) abort(glue::glue("Prior {type} must be strictly positive."))
  }
}

new_nma_prior <- function(dist, location = NA_real_, scale = NA_real_, df = NA_real_, ...) {
  o <- list(dist = dist,
            fun = deparse(sys.call(sys.parent())[[1]]),
            location = location,
            scale = scale,
            df = df,
            ...)
  class(o) <- "nma_prior"
  return(o)
}


#' Produce tidy prior details
#'
#' Produces prior details in a data frame, in a suitable format for
#' ggdist::stat_dist_*
#'
#' @param prior A nma_prior object
#' @param trunc Optional vector of truncation limits
#' @param ... Not used
#'
#' @return A data frame
#' @noRd
get_tidy_prior <- function(prior, trunc = NULL, ...) {
  if (!inherits(prior, "nma_prior"))
    abort("Not a `nma_prior` object.")
  if (!is.null(trunc) && (!rlang::is_double(trunc, n = 2) || trunc[2] <= trunc[1]))
    abort("`trunc` should be a numeric length 2 vector of lower and upper truncation limits.")

  d <- prior$dist
  is_trunc <- !is.null(trunc)

  if (d == "Normal") {
    out <- tibble::tibble(dist_label = d,
                          dist = "norm",
                          args = list(list(mean = prior$location, sd = prior$scale)))
  } else if (d == "half-Normal") {
    out <- tibble::tibble(dist_label = d,
                          dist = "trunc",
                          args = list(list(spec = "norm",
                                           a = if (is_trunc) max(0, trunc[1]) else 0,
                                           b = if (is_trunc) trunc[2] else Inf,
                                           mean = prior$location, sd = prior$scale)))
  } else if (d == "log-Normal") {
    out <- tibble::tibble(dist_label = d,
                          dist = "lnorm",
                          args = list(list(meanlog = prior$location, sdlog = prior$scale)))
  } else if (d == "Cauchy") {
    out <- tibble::tibble(dist_label = d,
                          dist = "cauchy",
                          args = list(list(location = prior$location, scale = prior$scale)))
  } else if (d == "half-Cauchy") {
    out <- tibble::tibble(dist_label = d,
                          dist = "trunc",
                          args = list(list(spec = "cauchy",
                                           a = if (is_trunc) max(0, trunc[1]) else 0,
                                           b = if (is_trunc) trunc[2] else Inf,
                                           location = prior$location, scale = prior$scale)))
  } else if (d == "Student t") {
    out <- tibble::tibble(dist_label = d,
                          dist = "gent",
                          args = list(list(location = prior$location, scale = prior$scale, df = prior$df)))
  } else if (d == "half-Student t") {
    out <- tibble::tibble(dist_label = d,
                          dist = "trunc",
                          args = list(list(spec = "gent",
                                           a = if (is_trunc) max(0, trunc[1]) else 0,
                                           b = if (is_trunc) trunc[2] else Inf,
                                           location = prior$location, scale = prior$scale, df = prior$df)))
  } else if (d == "log-Student t") {
    out <- tibble::tibble(dist_label = d,
                          dist = "logt",
                          args = list(list(location = prior$location, scale = prior$scale, df = prior$df)))
  } else if (d == "Exponential") {
    out <- tibble::tibble(dist_label = d,
                          dist = "exp",
                          args = list(list(rate = 1 / prior$scale)))
  } else if (d == "flat (implicit)") {
    out <- tibble::tibble(dist_label = d,
                          dist = "unif",
                          args = list(list(min = if (is_trunc) trunc[1] else -Inf,
                                           max = if (is_trunc) trunc[2] else Inf)))
  }

  if (is_trunc && d %in% c("Normal", "Cauchy", "Student t", "Exponential", "log-Normal", "log-Student t")) {
    out$args <- list(rlang::list2(spec = out$dist, a = trunc[1], b = trunc[2], !!! out$args[[1]]))
    out$dist <- "trunc"
  }

  return(out)
}


#' Get prior distribution call
#'
#' @param x A nma_prior object
#'
#' @return String giving call to construct x
#' @noRd
get_prior_call <- function(x) {
  # Deal with lists of priors
  if (all(purrr::map_lgl(x, ~inherits(., "nma_prior")))) {
    out <- purrr::imap_chr(x, ~paste0(.y, " = ", get_prior_call(.x)))
    out <- paste0("list(", paste(out, collapse = ", "), ")")
    return(out)
  }

  if (!inherits(x, "nma_prior"))
    abort("Not a `nma_prior` object.")

  prior_args <- purrr::list_modify(x, dist = purrr::zap(), fun = purrr::zap())
  prior_args <- prior_args[!is.na(prior_args)]
  if (length(prior_args)) {
    out <- paste0(x$fun, "(",
                  paste(names(prior_args), prior_args, sep = ' = ', collapse = ', '),
                  ")")
  } else {
    out <- paste0(x$fun, "()")
  }
  return(out)
}

# # Density for general truncated distribution
# dtrunc <- function(x, dist, trunc, ...) {
#   a <- list(...)
#   dfun <- paste0("d", dist)
#   pfun <- paste0("p", dist)
#
#   out <- ifelse(x < trunc[1] | x > trunc[2],
#                 0,
#                 do.call(dfun, args = rlang::list2(x = x, !!! a)) / (do.call(pfun, args = rlang::list2(q = trunc[2], !!! a)) - do.call(pfun, args = rlang::list2(q = trunc[1], !!! a))))
#
#   return(out)
# }
#
# # CDF for general truncated distribution
# ptrunc <- function(q, dist, trunc, ...) {
#   a <- list(...)
#   pfun <- paste0("p", dist)
#
#   out <- (do.call(pfun, args = rlang::list2(q = q, !!! a)) - do.call(pfun, args = rlang::list2(q = trunc[1], !!! a))) /
#     (do.call(pfun, args = rlang::list2(q = trunc[2], !!! a)) - do.call(pfun, args = rlang::list2(q = trunc[1], !!! a)))
#
#   return(out)
# }
#
# # Inverse CDF for general truncated distribution
# qtrunc <- function(p, dist, trunc, ...) {
#   a <- list(...)
#   qfun <- paste0("q", dist)
#   pfun <- paste0("p", dist)
#
#   pU <- do.call(pfun, args = rlang::list2(q = trunc[2], !!! a))
#   pL <- do.call(pfun, args = rlang::list2(q = trunc[1], !!! a))
#
#   pt <- pL + p * (pU - pL)
#
#   out <- do.call(qfun, args = rlang::list2(p = pt, !!! a))
#   out <- pmin(pmax(trunc[1], out), trunc[2])
#   return(out)
# }

#' Generalised Student's t distribution (with location and scale)
#'
#' Density, distribution, and quantile function for the generalised t
#' distribution with degrees of freedom `df`, shifted by `location` and scaled
#' by `scale`.
#'
#' @param x,q Vector of quantiles
#' @param p Vector of probabilities
#' @param df Degrees of freedom, greater than zero
#' @param location Location parameter
#' @param scale Scale parameter, greater than zero
#'
#' @return `dgent()` gives the density, `pgent()` gives the distribution
#'   function, `qgent()` gives the quantile function.
#' @export
#' @rdname generalised_t
#'
dgent <- function(x, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")

  return(stats::dt((x - location) / scale, df) / scale)
}

#' @export
#' @rdname generalised_t
pgent <- function(q, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")

  return(stats::pt((q - location) / scale, df))
}

#' @export
#' @rdname generalised_t
qgent <- function(p, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")

  return(location + scale * stats::qt(p, df))
}

#' Log Student's t distribution
#'
#' Density, distribution, and quantile function for the log t distribution,
#' whose logarithm has degrees of freedom `df`, mean `location`, and standard
#' deviation `scale`.
#'
#' @param x,q Vector of quantiles
#' @param p Vector of probabilities
#' @param df Degrees of freedom, greater than zero
#' @param location Location parameter
#' @param scale Scale parameter, greater than zero
#'
#' @details If \eqn{\log(Y) \sim t_\nu(\mu, \sigma^2)}, then \eqn{Y} has a log t
#'   distribution with `location` \eqn{\mu}, `scale` \eqn{\sigma}, and `df`
#'   \eqn{\nu}.
#'
#'   The mean and all higher moments of the log t distribution are undefined or
#'   infinite.
#'
#'   If `df = 1` then the distribution is a log Cauchy distribution. As `df`
#'   tends to infinity, this approaches a log Normal distribution.
#'
#' @return `dlogt()` gives the density, `plogt()` gives the distribution
#'   function, `qlogt()` gives the quantile function.
#' @export
#' @rdname log_t
#'
dlogt <- function(x, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")
  if (!is.numeric(df) || any(df <= 0))
    abort("`df` must be numeric, greater than zero.")

  return(dgent(log(x), df = df, location = location, scale = scale) / x)
}

#' @export
#' @rdname log_t
plogt <- function(q, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")

  return(pgent(q = log(q), df = df, location = location, scale = scale))
}

#' @export
#' @rdname log_t
qlogt <- function(p, df, location = 0, scale = 1) {
  if (!is.numeric(location))
    abort("`location` must be numeric.")
  if (!is.numeric(scale) || any(scale <= 0))
    abort("`scale` must be numeric, greater than zero.")

  return(exp(qgent(p = p, df = df, location = location, scale = scale)))
}
