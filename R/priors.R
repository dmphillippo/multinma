#' Prior distributions
#'
#' These functions are used to specify prior distributions for the model parameters.
#'
#' @param mean Prior mean.
#' @param median Prior median.
#' @param sd Prior standard deviation.
#' @param scale Prior scale.
#' @param df Prior degrees of freedom.
#' @param rate Prior rate.
#' @param lower,upper Lower and upper bounds for a uniform prior distribution.
#'
#' @rdname priors
#' @aliases normal
#'
#' @return Object of class [prior].
#' @export
#'
#' @examples
normal <- function(mean = 0, sd) {
  check_prior_location(mean)
  check_prior_scale(sd)
}

#' @rdname priors
#' @export
half_normal <- function(sd) {
  check_prior_scale(sd)
}

#' @rdname priors
#' @export
cauchy <- function(median = 0, scale) {
  check_prior_location(median, type = "median")
  check_prior_scale(scale, type = "scale")
}

#' @rdname priors
#' @export
half_cauchy <- function(scale) {
  check_prior_scale(scale, type = "scale")
}

#' @rdname priors
#' @export
student_t <- function(mean = 0, sd, df) {
  check_prior_location(mean)
  check_prior_scale(sd)
  check_prior_scale(df, type = "degrees of freedom")
}

#' @rdname priors
#' @export
half_student_t <- function(sd, df) {
  check_prior_scale(sd)
  check_prior_scale(df, type = "degrees of freedom")
}

#' @rdname priors
#' @export
exponential <- function(rate) {
  check_prior_scale(rate, type = "rate")
}

#' @rdname priors
#' @export
uniform <- function(lower, upper) {
  check_prior_location(lower, type = "lower limit")
  check_prior_location(lower, type = "upper limit")
}

# Check functions
check_prior_location <- function(x, type = "mean") {
  if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
  if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
}

check_prior_scale <- function(x, type = "standard deviation") {
  if (!is.numeric(x)) abort(glue::glue("Prior {type} must be numeric."))
  if (length(x) > 1) abort(glue::glue("Prior {type} must be numeric, length 1."))
  if (x <= 0) abort(glue::glue("Prior {type} must be strictly positive."))
}
