#' Bayesian R-squared
#'
#' Calculate Bayesian-$R^2$ values for a regression model. Currently this is
#' only calculated for IPD regression models. Any models with both IPD and AgD
#' (e.g. ML-NMR) will ignore the AgD part for this calculation.
#'
#' @param object A `stan_nma` object.
#' @param ... Not used.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a 3D MCMC array
#'  of samples of the Bayesian $R^2$.
#' @aliases bayes_R2
#' @importFrom rstantools bayes_R2
#' @export
bayes_R2.stan_nma <- function(object, ..., probs = c(0.025, 0.5, 0.975), summary = TRUE) {

  if (!rlang::is_bool(summary)) abort("`summary` must be TRUE or FALSE.")

  if (!has_ipd(object$network)) abort("No IPD present for which to calculate R2.")
  if (has_agd_arm(object$network) || has_agd_contrast(object$network))
    inform("Note: R-squared calculated on IPD portion of the model only.")
  if (object$likelihood %in% valid_lhood$survival) abort("Not supported for survival outcomes.")

  # Array of predicted responses
  mu_pred <- as.array(object, pars = "fitted_ipd")

  odim <- c(dim(mu_pred)[1:2], 1)
  onames <- list(iterations = NULL, chains = NULL, parameters = "bayes_R2")
  var_res <- var_mu_pred <- array(NA, dim = odim, dimnames = onames)

  # Model variance of residuals
  if (object$likelihood %in% c(valid_lhood$binary, valid_lhood$count)) {
    var_res[,,1] <- apply(mu_pred * (1 - mu_pred), 1:2, "mean")
  } else if (object$likelihood == "ordered") {
    n <- nrow(object$network$ipd)
    var_res[,,1] <- apply(mu_pred, 1:2, function(p) {
      pmat <- matrix(p, nrow = n, byrow = TRUE)[, -1, drop = FALSE]
      mean(rowSums(pmat * (1 - pmat)))
    })
  } else if (object$likelihood == "normal") {
    ss <- dplyr::group_by(object$network$ipd, .data$.study, .data$.trt) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::arrange(.data$.study, .data$.trt) %>%
      dplyr::pull("n")
    var_res[,,1] <- apply(as.array(object, pars = "sigma")^2, 1:2, FUN = weighted.mean, w = ss)
  } else if (object$likelihood == "poisson") {
    var_res[,,1] <- apply(mu_pred, 1:2, "mean")
  } else {
    abort(glue::glue("Likelihood '{object$likelihood}' not yet supported."))
  }

  if (object$likelihood == "ordered") {
    var_mu_pred[,,1] <- apply(mu_pred, 1:2, function(p) {
      pmat <- matrix(p, nrow = n, byrow = TRUE)[, -1, drop = FALSE]
      sum(diag(var(pmat)))
    })
  } else {
    var_mu_pred[,,1] <- apply(mu_pred, 1:2, "var")
  }


  r2 <- var_mu_pred / (var_mu_pred + var_res)

  if (summary) {
    out <- list(summary = summary_mcmc_array(r2, probs), sims = r2)
    class(out) <- "nma_summary"
  } else {
    out <- r2
  }

  return(out)
}

#' @rdname bayes_R2.stan_nma
#' @aliases loo_R2
#' @importFrom rstantools loo_R2
#' @export
loo_R2.stan_nma <- function(object, ..., probs = c(0.025, 0.5, 0.975), summary = TRUE) {

  require_pkg("loo")

  if (!rlang::is_bool(summary)) abort("`summary` must be TRUE or FALSE.")

  if (!has_ipd(object$network)) abort("No IPD present for which to calculate R2.")
  if (has_agd_arm(object$network) || has_agd_contrast(object$network))
    inform("Note: R-squared calculated on IPD portion of the model only.")
  if (object$likelihood %in% valid_lhood$survival) abort("Not supported for survival outcomes.")

  if (object$likelihood == "ordered") abort("Not supported for ordered outcomes.")

  # Observed outcomes (non-ordered, unchanged)
  if (object$likelihood %in% c(valid_lhood$binary, valid_lhood$count, "poisson")) {
    y <- object$network$ipd$.r
  } else if (object$likelihood == "normal") {
    y <- object$network$ipd$.y
  } else {
    abort(glue::glue("Likelihood '{object$likelihood}' not yet supported."))
  }

  n <- length(y)

  # Array of predicted responses
  mu_pred <- as.array(object, pars = "fitted_ipd")

  # PSIS (matrix form required from here)
  mu_pred_mat <- as.matrix.nma_summary(mu_pred)

  log_ratios <- -as.matrix(object, pars = "log_lik")[, 1:n, drop = FALSE]
  psis_object <- loo::psis(log_ratios)

  mu_pred_loo <- loo::E_loo(mu_pred_mat, psis_object, log_ratios = log_ratios)$value
  err_loo <- mu_pred_loo - y

  # Dirichlet weights Bayesian bootstrap
  niter <- nrow(mu_pred_mat)
  exp_draws <- matrix(rexp(n * niter, rate = 1), nrow = niter, ncol = n)
  wts <- exp_draws / rowSums(exp_draws)

  var_y <- (rowSums(sweep(wts, 2, y^2, FUN = "*")) -
              rowSums(sweep(wts, 2, y, FUN = "*"))^2) * (n/(n-1))

  var_err_loo <- (rowSums(sweep(wts, 2, err_loo^2, FUN = "*")) -
                    rowSums(sweep(wts, 2, err_loo, FUN = "*")^2)) * (n/(n-1))


  odim <- c(dim(mu_pred)[1:2], 1)
  onames <- list(iterations = NULL, chains = NULL, parameters = "loo_R2")
  r2 <- array(NA, dim = odim, dimnames = onames)

  r2[,,1] <- 1 - var_err_loo / var_y
  r2[r2 < -1] <- -1
  r2[r2 > 1] <- 1

  if (summary) {
    out <- list(summary = summary_mcmc_array(r2, probs), sims = r2)
    class(out) <- "nma_summary"
  } else {
    out <- r2
  }

  return(out)
}

