#' Prior distributions
#'
#' These functions are used to specify prior distributions for the model parameters.
#'
#' @param location Prior location. Typically prior mean, except for Cauchy which is prior median.
#' @param scale Prior scale. Typically prior standard deviation.
#' @param df Prior degrees of freedom.
#' @param rate Prior rate.
# #' @param lower,upper Lower and upper bounds for a uniform prior distribution.
#'
#' @rdname priors
#' @name priors
#' @aliases priors normal
#'
#' @return Object of class [nma_prior].
#' @export
#'
#' @examples
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
exponential <- function(scale = 1/rate, rate = 1/scale) {
  check_prior_scale(rate, type = "scale or rate")

  return(new_nma_prior("Exponential", scale = scale))
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
  if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
  if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
}

check_prior_scale <- function(x, type = "scale (standard deviation)") {
  if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
  if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
  if (x <= 0) abort(glue::glue("Prior {type} must be strictly positive."))
}

new_nma_prior <- function(dist, location = NA_real_, scale = NA_real_, df = NA_real_, ...) {
  o <- list(dist = dist,
            location = location,
            scale = scale,
            df = df,
            ...)
  class(o) <- "nma_prior"
  return(o)
}
