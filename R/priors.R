#' Prior distributions
#'
#' These functions are used to specify prior distributions for the model parameters.
#'
#' @param mean Prior mean.
#' @param median Prior median.
#' @param sd Prior standard deviation.
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

}

#' @rdname priors
#' @export
half_normal <- function(sd) {

}

#' @rdname priors
#' @export
cauchy <- function(median = 0, sd) {

}

#' @rdname priors
#' @export
student_t <- function(mean = 0, sd, df) {

}

#' @rdname priors
#' @export
exponential <- function(rate) {

}

#' @rdname priors
#' @export
uniform <- function(lower, upper) {

}

# Check functions
check_prior_location <- function(x, type = "mean") {

}

check_prior_scale <- function(x, type = "sd") {

}
