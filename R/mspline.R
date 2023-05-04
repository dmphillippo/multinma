#' Distribution functions for M-spline baseline hazards
#'
#' Density, distribution, quantile, hazard, cumulative hazard, restricted mean
#' survival time and mean survival time functions for the M-spline baseline
#' hazards model.
#'
#' @param x,q Vector of quantiles
#' @param basis M-spline basis produced by [splines2::mSpline()] or
#'   [splines2::iSpline()]
#' @param scoef Vector (or matrix) of spline coefficients with length (or number
#'   of columns) equal to the dimension of `basis`
#' @param rate Vector of rate parameters
#' @param log,log.p Logical; if `TRUE`, probabilities `p` are given as
#'   \eqn{\log(p)}
#' @param lower.tail Logical; if `TRUE` (the default), probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}
#' @param p Vector of probabilities
#' @param t Vector of times to which the restricted mean survival time is
#'   calculated
#' @param start Optional left-truncation time or times. The returned restricted
#'   mean survival will be conditioned on survival up to this time
#'
#' @details Survival models with a flexible M-spline on the baseline hazard are
#'   described by \insertRef{Brilleman2020}{multinma}. Piecewise-exponential
#'   baseline hazards are a special case where the degree of the M-spline
#'   polynomial is 0.
#'
#'   The d/p/h/H functions are calculated from their definitions. `qmspline()`
#'   uses numerical inversion via [flexsurv::qgeneric()]. `rmst_mspline()` and
#'   `mean_mspline()` use numerical integration via [flexsurv::rmst_generic()].
#'
#' @return `dmspline()` gives the density, `pmspline()` gives the distribution
#'   function (CDF), `qmspline()` gives the quantile function (inverse-CDF),
#'   `hmspline()` gives the hazard function, `Hmspline()` gives the cumulative
#'   hazard function, `rmst_mspline()` gives restricted mean survival times, and
#'   `mean_mspline()` gives mean survival times.
#'
#' @rdname mspline
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
dmspline <- function(x, basis, scoef, rate, log = FALSE) {
  if (!inherits(basis, "mSpline")) abort("`basis` must be an M-spline basis produced by splines2::mSpline() or splines2::iSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  if (missing(x)) out <- hmspline(basis = basis, scoef = scoef, rate = rate) * pmspline(basis = basis, scoef = scoef, rate = rate, lower.tail = FALSE)
  else out <- hmspline(x, basis = basis, scoef = scoef, rate = rate) * pmspline(q = x, basis = basis, scoef = scoef, rate = rate, lower.tail = FALSE)

  if (log) out <- log(out)

  return(out)
}

#' @rdname mspline
#' @export
pmspline <- function(q, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE) {
  if (!inherits(basis, "mSpline")) abort("`basis` must be an M-spline basis produced by splines2::mSpline() or splines2::iSpline()")
  if (!rlang::is_bool(lower.tail)) abort("`lower.tail` must be a logical value (TRUE or FALSE).")
  if (!rlang::is_bool(log.p)) abort("`log.p` must be a logical value (TRUE or FALSE).")

  if (missing(q)) out <- exp(-Hmspline(basis = basis, scoef = scoef, rate = rate))
  else  out <- exp(-Hmspline(x = q, basis = basis, scoef = scoef, rate = rate))

  if (lower.tail) out <- 1 - out
  if (log.p) out <- log(out)

  return(out)
}

#' @rdname mspline
#' @export
qmspline <- function(p, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE) {
  if (!is.numeric(p)) abort("`p` must be a numeric vector of quantiles.")

  flexsurv::qgeneric(pmspline, p = p,
                     scalarargs = "basis", matargs = "scoef",
                     basis = basis, scoef = scoef, rate = rate,
                     lower.tail = lower.tail, log.p = log.p,
                     lbound = 0)
}

#' @rdname mspline
#' @export
hmspline <- function(x, basis, scoef, rate, log = FALSE) {
  if (!inherits(basis, "mSpline")) abort("`basis` must be an M-spline basis produced by splines2::mSpline() or splines2::iSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  if (is.matrix(scoef)) {
    if (missing(x)) out <- mapply(crossprod, asplit(update(basis, integral = FALSE), 1), asplit(scoef, 1)) * rate
    else  out <- mapply(crossprod, asplit(update(basis, x = x, integral = FALSE), 1), asplit(scoef, 1)) * rate
  } else {
    if (missing(x)) out <- update(basis, integral = FALSE) %*% scoef * rate
    else  out <- update(basis, x = x, integral = FALSE) %*% scoef * rate
  }

  # Return 0 for x < 0
  lt0 <- x < 0
  out[lt0] <- 0

  if (log) out <- log(out)

  return(drop(out))
}

#' @rdname mspline
#' @export
Hmspline <- function(x, basis, scoef, rate, log = FALSE) {
  if (!inherits(basis, "mSpline")) abort("`basis` must be an M-spline basis produced by splines2::mSpline() or splines2::iSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  if (is.matrix(scoef)) {
    if (missing(x)) out <- mapply(crossprod, asplit(update(basis, integral = TRUE), 1), asplit(scoef, 1)) * rate
    else  out <- mapply(crossprod, asplit(update(basis, x = x, integral = TRUE), 1), asplit(scoef, 1)) * rate
  } else {
    if (missing(x)) out <- update(basis, integral = TRUE) %*% scoef * rate
    else  out <- update(basis, x = x, integral = TRUE) %*% scoef * rate
  }

  # Return 0 for x < 0
  lt0 <- x < 0
  out[lt0] <- 0

  if (log) out <- log(out)

  return(drop(out))
}

#' @rdname mspline
#' @export
rmst_mspline <- function(t, basis, scoef, rate, start = 0) {
  flexsurv::rmst_generic(pmspline, t, start = start,
                         basis = basis, scoef = scoef, rate = rate,
                         scalarargs = "basis", matargs = "scoef")
}

#' @rdname mspline
#' @export
mean_mspline <- function(basis, scoef, rate, ...) {
  rmst_mspline(t = Inf, basis, scoef, rate)
}
