#' Distribution functions for M-spline baseline hazards
#'
#' Density, distribution, quantile, hazard, cumulative hazard, and restricted
#' mean survival time functions for the M-spline baseline hazards model.
#'
#' @param x,q Vector of quantiles
#' @param basis M-spline basis produced by [splines2::mSpline()]
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
#'   described by \insertCite{Brilleman2020;textual}{multinma}.
#'   Piecewise-exponential baseline hazards are a special case where the degree
#'   of the M-spline polynomial is 0.
#'
#'   The d/p/h/H functions are calculated from their definitions. `qmspline()`
#'   uses numerical inversion via [flexsurv::qgeneric()]. `rmst_mspline()`uses
#'   numerical integration via [flexsurv::rmst_generic()], except for the
#'   special case of the piecewise-exponential hazard (i.e. degree 0 M-splines)
#'   which uses the explicit formula from
#'   \insertCite{Royston2013;textual}{multinma}.
#'
#'   Beyond the boundary knots, the hazard is assumed to be constant. (This
#'   differs from the approach in [splines2::mSpline()] that extrapolates the
#'   polynomial basis functions, which is numerically unstable and highly
#'   dependent on the data just before the boundary knots.) As with all
#'   extrapolation, care should be taken when evaluating the splines at times
#'   beyond the boundary knots (either directly through the d/p/h/H/rmst
#'   functions, or indirectly by requesting quantiles with `qmspline()` that
#'   correspond to times beyond the boundary knots). For this reason evaluating
#'   the (unrestricted) mean survival time is not generally recommended as this
#'   requires integrating over an infinite time horizon (i.e. `rmst_mspline()`
#'   with `t = Inf`).
#'
#' @return `dmspline()` gives the density, `pmspline()` gives the distribution
#'   function (CDF), `qmspline()` gives the quantile function (inverse-CDF),
#'   `hmspline()` gives the hazard function, `Hmspline()` gives the cumulative
#'   hazard function, and `rmst_mspline()` gives restricted mean survival times.
#'
#' @rdname mspline
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
dmspline <- function(x, basis, scoef, rate, log = FALSE) {
  if (!is_mspline(basis)) abort("`basis` must be an M-spline basis produced by splines2::mSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  if (missing(x)) out <- hmspline(basis = basis, scoef = scoef, rate = rate) * pmspline(basis = basis, scoef = scoef, rate = rate, lower.tail = FALSE)
  else out <- hmspline(x, basis = basis, scoef = scoef, rate = rate) * pmspline(q = x, basis = basis, scoef = scoef, rate = rate, lower.tail = FALSE)

  if (log) out <- log(out)

  return(out)
}

#' @rdname mspline
#' @export
pmspline <- function(q, basis, scoef, rate, lower.tail = TRUE, log.p = FALSE) {
  if (!is_mspline(basis)) abort("`basis` must be an M-spline basis produced by splines2::mSpline()")
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
  require_pkg("flexsurv")

  flexsurv::qgeneric(pmspline, p = p,
                     scalarargs = "basis", matargs = "scoef",
                     basis = basis, scoef = scoef, rate = rate,
                     lower.tail = lower.tail, log.p = log.p,
                     lbound = 0)
}

#' @rdname mspline
#' @export
hmspline <- function(x, basis, scoef, rate, log = FALSE) {
  require_pkg("splines2")
  if (!is_mspline(basis)) abort("`basis` must be an M-spline basis produced by splines2::mSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  # Extrapolate with constant hazard beyond boundary knots
  if (!missing(x)) {
    x <- pmax(x, attr(basis, "Boundary.knots")[1])
    x <- pmin(x, attr(basis, "Boundary.knots")[2])
  }

  xb <- if (missing(x)) update(basis, integral = FALSE)
  else update(basis, x = x, integral = FALSE)

  if (!is.matrix(scoef)) scoef <- matrix(scoef, nrow = 1)

  if (is.matrix(scoef) && nrow(scoef) == nrow(xb)) {
    out <- rowSums(xb * scoef) * rate
  } else if (nrow(xb) == 1) {
    out <- scoef %*% xb[1,] * rate
  } else if (is.matrix(scoef) && nrow(scoef) == 1) {
    out <- xb %*% scoef[1,] * rate
  } else {
    out <- xb %*% scoef * rate
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
  require_pkg("splines2")
  if (!is_mspline(basis)) abort("`basis` must be an M-spline basis produced by splines2::mSpline()")
  if (!rlang::is_bool(log)) abort("`log` must be a logical value (TRUE or FALSE).")

  # Extrapolate with constant hazard beyond boundary knots
  if (!missing(x)) {
    lower <- attr(basis, "Boundary.knots")[1]
    upper <- attr(basis, "Boundary.knots")[2]

    ex_lo <- x < lower & x > 0
    ex_up <- x > upper

    ex_time_lo <- x[ex_lo]
    ex_time_up <- x[ex_up] - upper

    xorig <- x
    x <- pmax(x, lower)
    x <- pmin(x, upper)
  } else {
    ex_lo <- ex_up <- FALSE
  }

  xb <- if (missing(x)) update(basis, integral = TRUE)
        else update(basis, x = x, integral = TRUE)

  if (!is.matrix(scoef)) scoef <- matrix(scoef, nrow = 1)

  if (is.matrix(scoef) && nrow(scoef) == nrow(xb)) {
    out <- rowSums(xb * scoef) * rate
  } else if (nrow(xb) == 1) {
    out <- scoef %*% xb[1,] * rate
  } else if (is.matrix(scoef) && nrow(scoef) == 1) {
    out <- xb %*% scoef[1,] * rate
  } else {
    out <- xb %*% scoef * rate
  }

  if (any(ex_up)) {
    ex_up_rate <- if (length(rate) == 1) rate else rate[ex_up]
    ex_up_scoef <- if (nrow(scoef) == 1) scoef else scoef[ex_up, , drop = FALSE]
    h_up <- hmspline(xorig[ex_up], basis = basis, scoef = ex_up_scoef, rate = ex_up_rate)
    out[ex_up] <- out[ex_up] + h_up * ex_time_up
  }

  if (any(ex_lo)) {
    ex_lo_rate <- if (length(rate) == 1) rate else rate[ex_lo]
    ex_lo_scoef <- if (nrow(scoef) == 1) scoef else scoef[ex_lo, , drop = FALSE]
    h_lo <- hmspline(xorig[ex_lo], basis = basis, scoef = ex_lo_scoef, rate = ex_lo_rate)
    out[ex_lo] <- h_lo * ex_time_lo
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
  if (attr(basis, "degree") == 0) {
    # Piecewise exponential has explicit formula (Royston and Parmar 2013)
    nr <- max(length(t), if (is.matrix(scoef)) nrow(scoef) else 1, length(rate))

    knots <- c(attr(basis, "Boundary.knots")[1], attr(basis, "knots"))
    knots <- matrix(knots, nrow = nr, ncol = length(knots), byrow = TRUE)

    h <- apply(knots, 2, hmspline, basis = basis, scoef = scoef, rate = rate)

    delta <- t(apply(pmax(pmin(cbind(knots, Inf), t) - start, 0), 1, diff))

    hd <- h * delta

    H <- t(apply(cbind(0, hd[, -ncol(knots), drop = FALSE]), 1, cumsum))

    rowSums(exp(-H) / h * (1 - exp(-hd)))

  } else {
    # General M-splines require numerical integration
    require_pkg("flexsurv")
    flexsurv::rmst_generic(pmspline, t, start = start,
                           basis = basis, scoef = scoef, rate = rate,
                           scalarargs = "basis", matargs = "scoef")
  }
}

# Don't export mean_mspline - this is a bad idea in general
#' @noRd
mean_mspline <- function(basis, scoef, rate, ...) {
  rmst_mspline(t = Inf, basis, scoef, rate)
}

# Check for mspline/ispline objects
#' @param x Object to check is MSpline class
#' @return logical
#' @noRd
is_mspline <- function(x) {
  inherits(x, "MSpline")
}

# Create a (logit scale) coefficient vector corresponding to a constant hazard
#' @param basis A M-spline basis created using splines2::mSpline()
#' @return Coefficient vector
#' @noRd
mspline_constant_hazard <- function(basis) {
  if (!is_mspline(basis)) abort("`basis` must be an M-spline basis created using splines2::mSpline()")

  df <- ncol(basis)
  ord <- attr(basis, "degree") + 1
  iknots <- attr(basis, "knots")
  bknots <- attr(basis, "Boundary.knots")

  # Using approach of Jackson arXiv:2306.03957
  knots <- c(rep(bknots[1], ord), iknots, rep(bknots[2], ord))
  coefs <- (knots[(1:df) + ord] - knots[1:df]) / (ord * (diff(bknots)))

  # inverse softmax transform
  inv_softmax(coefs)
}

#' Create a vector of (normalised) weights for a RW(1) prior on an M-spline with
#' unequally spaced knots
#' @param basis M-spline basis created using splines2::mSpline()
#' @return Vector of normalised weights
#' @noRd
rw1_prior_weights <- function(basis) {
  nscoef <- ncol(basis)
  ord <- attr(basis, "degree") + 1
  iknots <- attr(basis, "knots")
  bknots <- attr(basis, "Boundary.knots")
  knots <- c(rep(bknots[1], ord), iknots, rep(bknots[2], ord))
  if (ord == 1) {
    wts <- (knots[2:nscoef] - knots[1:(nscoef - 1)]) / (knots[nscoef] - bknots[1])
  } else {
    wts <- (knots[(ord + 1):(nscoef + ord - 1)] - knots[2:nscoef]) / ((ord - 1) * (bknots[2] - bknots[1]))
  }
  return(sqrt(wts))
}

#' Softmax transform
#'
#' The softmax transform is a multivariate generalisation of the logit
#' transform. `softmax()` maps a vector of \eqn{K-1} values on the real line to a
#' \eqn{K}-simplex (i.e. values between 0 and 1, that sum to 1). `inv_softmax()`
#' provides the inverse transform, mapping a \eqn{K}-simplex vector to a vector of
#' \eqn{K-1} real values.
#'
#' @param x \eqn{K-1} vector of reals
#' @return `softmax()` returns a vector of length \eqn{K} that is a simplex.
#'   `inv_softmax()` returns a vector of reals of length \eqn{K-1}.
#' @export
#'
#' @examples
#' x <- c(-1, 3, -0.5, 2)
#' (p <- softmax(x))
#' sum(p)
#' inv_softmax(p)
#'
softmax <- function(x) {
  x0 <- c(0, x)
  exp(x0 - logsumexp(x0))
}

logsumexp <- function(x) {
  maxx <- max(x)
  max(x) + log(sum(exp(x - maxx)))
}

#' inverse softmax transform
#' @param p \eqn{K} vector simplex
#' @export
#' @rdname softmax
inv_softmax <- function(p) {
  log(p[-1]) - log(p[1])
}


#' Knot locations for M-spline baseline hazard models
#'
#' Several different algorithms are provided to calculate knot locations for
#' M-spline baseline hazard models. This function is called internally within
#' the [nma()] function, but may be called directly by the user for more
#' control.
#'
#' @param network A network object, containing survival outcomes
#' @param n_knots Non-negative integer giving the number of internal knots
#'   (default `7`)
# #' @param degree Non-negative integer giving the degree of the M-spline
# #'   polynomial (default `3`)
#' @param type String specifying the knot location algorithm to use (see
#'   details). The default used by [nma()] is `"quantile"`, except when a
#'   regression model is specified (using `aux_regression`) in which case the
#'   default is `"quantile_common"`.
#'
#' @details The `type` argument can be used to choose between different
#'   algorithms for placing the knots:
#'
#' \describe{
#'   \item{`"quantile"`}{Creates separate knot locations for each study,
#'   internal knots are placed at evenly-spaced quantiles of the observed event
#'   times within each study.}
#'   \item{`"quantile_lumped"`}{Creates a common set of knots for all studies,
#'   calculated as evenly-spaced quantiles of the observed event times from all
#'   studies lumped together.}
#'   \item{`"quantile_common"`}{Creates a common set of knots for all studies,
#'   taking quantiles of the quantiles of the observed event times within each
#'   study. This often seems to result in a more even knot spacing than
#'   `"quantile_lumped"`, particularly when follow-up is uneven across studies,
#'   and may handle differing behaviour in the baseline hazard across studies
#'   better than `"quantile_longest"`.}
#'   \item{`"quantile_longest"`}{Creates a common set of knots for all studies,
#'   using evenly-spaced quantiles of the observed event times in the longest
#'   study.}
#'   \item{`"equal"`}{Creates separate knot locations for each study, at
#'   evenly-spaced times between the boundary knots in each study.}
#'   \item{`"equal_common"`}{Creates a common set of knots for all studies, at
#'   evenly-spaced times between the earliest entry time and last
#'   event/censoring time in the network.}
#' }
#'
#' Boundary knots are calculated as follows:
#'
#' * For separate knot locations in each study, boundary knots are placed at the
#' earliest entry time and last event/censoring time in each study.
#' * For a common set of knots across all studies, boundary knots are placed at
#' the earliest entry time and last event/censoring time across all studies.
#'
#' Models with regression on the spline coefficients (i.e. with `aux_regression`
#' specified) require a common set of knots across all studies.
#'
#' Provided that a sufficient number of knots are used, model fit should be
#' largely unaffected by the knot locations. However, sampling difficulties can
#' sometimes occur if knot placement is poor, for example if a knot is placed
#' just before the last follow-up time in a study.
#'
#' @return A named list of vectors giving the knot locations in each study.
#' @export
#'
#' @template ex_ndmm_network
#' @examples
#'
#' # The default knot locations
#' make_knots(ndmm_net, type = "quantile")
#'
#' # Increasing the number of knots
#' make_knots(ndmm_net, n_knots = 10)
#'
#' # Comparing alternative knot positioning algorithms
#' # Visualise these with a quick function
#' plot_knots <- function(network, knots) {
#'   ggplot2::ggplot() +
#'     geom_km(network) +
#'     ggplot2::geom_vline(ggplot2::aes(xintercept = .data$knot),
#'                         data = tidyr::pivot_longer(as.data.frame(knots), cols = dplyr::everything(),
#'                                                    names_to = "Study", values_to = "knot"),
#'                         linetype = 2, colour = "grey60") +
#'     ggplot2::facet_wrap(~Study) +
#'     theme_multinma()
#' }
#'
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "quantile"))
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "quantile_common"))
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "quantile_lumped"))
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "quantile_longest"))
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "equal"))
#' plot_knots(ndmm_net, make_knots(ndmm_net, type = "equal_common"))
#'
make_knots <- function(network,
                       n_knots = 7,
                       type = c("quantile",
                                "quantile_common",
                                "quantile_lumped",
                                "quantile_longest",
                                "equal",
                                "equal_common")) {

  # Argument checks
  if (!inherits(network, "nma_data"))
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")

  if (!(identical(network$outcome$ipd, "survival") || identical(network$outcome$agd_arm, "survival")))
    abort("`network` does not contain survival outcomes.")

  if (!rlang::is_scalar_integerish(n_knots, finite = TRUE) || n_knots < 1)
    abort("`n_knots` must be a positive integer giving the number of internal knots.")

  type <- rlang::arg_match(type)

  # Get survival data
  survdat <- dplyr::bind_rows(
    if (has_ipd(network)) dplyr::select(network$ipd, ".study", ".trt", ".Surv") else NULL,
    if (has_agd_arm(network)) dplyr::select(network$agd_arm, ".study", ".trt", ".Surv") %>% tidyr::unnest(cols = ".Surv") else NULL
  )

  survdat <- dplyr::mutate(survdat,
                           !!! get_Surv_data(survdat$.Surv),
                           observed = .data$status == 1)

  observed_survdat <- dplyr::filter(survdat, .data$observed)

  studies <- unique(survdat$.study)
  n_studies <- length(studies)

  # Calculate knots
  if (type == "quantile") {

    # Boundary knots
    b_knots <- by(survdat, survdat$.study, function(x) c(min(x$delay_time), max(x$time)),
                  simplify = FALSE)

    # Internal knots
    i_knots <- by(observed_survdat, observed_survdat$.study,
                  function(x) quantile(x$time, probs = seq(1, n_knots) / (n_knots + 1)),
                  simplify = FALSE)

  } else if (type == "quantile_common") {

    # Boundary knots
    b_knots <- c(min(survdat$delay_time), max(survdat$time))
    b_knots <- rep_len(list(b_knots), n_studies)
    names(b_knots) <- studies

    # Internal knots
    # Take quantiles of knots and boundary knots in individual studies
    temp_b_knots <- by(survdat, survdat$.study, function(x) c(min(x$delay_time), max(x$time)),
                       simplify = FALSE)

    temp_i_knots <- by(observed_survdat, observed_survdat$.study,
                       function(x) quantile(x$time, probs = seq(1, n_knots) / (n_knots + 1)),
                       simplify = FALSE)

    i_knots <- quantile(c(unlist(temp_i_knots), unlist(temp_b_knots)), probs = seq(0, 1, length.out = n_knots+2))[2:(n_knots+1)]
    i_knots <- rep_len(list(i_knots), n_studies)
    names(i_knots) <- studies

  } else if (type == "quantile_lumped") {

    # Boundary knots
    b_knots <- c(min(survdat$delay_time), max(survdat$time))
    b_knots <- rep_len(list(b_knots), n_studies)
    names(b_knots) <- studies

    # Internal knots
    i_knots <- quantile(observed_survdat$time, probs = seq(1, n_knots) / (n_knots + 1))
    i_knots <- rep_len(list(i_knots), n_studies)
    names(i_knots) <- studies

  } else if (type == "quantile_longest") {

    # Boundary knots
    b_knots <- c(min(survdat$delay_time), max(survdat$time))
    b_knots <- rep_len(list(b_knots), n_studies)
    names(b_knots) <- studies

    # Find study with longest follow-up
    fu <- by(survdat, survdat$.study, function(x) max(x$time) - min(x$delay_time), simplify = TRUE)
    max_fu <- names(which.max(fu))

    # Internal knots
    i_knots <- quantile(dplyr::filter(observed_survdat, .data$.study == max_fu)$time, probs = seq(1, n_knots) / (n_knots + 1))
    i_knots <- rep_len(list(i_knots), n_studies)
    names(i_knots) <- studies

  } else if (type == "equal") {

    # Boundary knots
    b_knots <- by(survdat, survdat$.study, function(x) c(min(x$delay_time), max(x$time)),
                  simplify = FALSE)

    # Internal knots
    i_knots <- by(survdat, survdat$.study,
                  function(x) seq(from = min(x$delay_time), to = max(x$time), length.out = n_knots+2)[2:(n_knots+1)],
                  simplify = FALSE)

  } else if (type == "equal_common") {

    # Boundary knots
    b_knots <- c(min(survdat$delay_time), max(survdat$time))
    b_knots <- rep_len(list(b_knots), n_studies)
    names(b_knots) <- studies

    # Internal knots
    i_knots <- seq(from = min(survdat$delay_time), to = max(survdat$time), length.out = n_knots+2)[2:(n_knots+1)]
    i_knots <- rep_len(list(i_knots), n_studies)
    names(i_knots) <- studies
  }

  # Combine boundary and internal knots
  out <- purrr::map2(b_knots, i_knots, ~ c(.x[1], .y, .x[2]))

  return(out)
}
