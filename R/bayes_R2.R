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
    # This is still weird...
    n <- nrow(object$network$ipd)
    var_res[,,1] <- apply(mu_pred, 1:2, function(p) mean(apply(matrix(p, nrow = n), 1, "prod")))
    # var_res[,,1] <- apply(mu_pred, 1:2, function(p) mean(apply(matrix(p*(1-p), nrow = n)[, -1, drop = FALSE], 1, "sum")))
  } else if (object$likelihood == "normal") {
    var_res[,,1] <- as.array(object, pars = "sigma")^2
  } else if (object$likelihood == "poisson") {
    var_res[,,1] <- apply(mu_pred, 1:2, "mean")
  } else {
    abort(glue::glue("Likelihood '{object$likelihood}' not yet supported."))
  }

  if (object$likelihood == "ordered") {
    # This is still weird...
    var_mu_pred[,,1] <- apply(mu_pred, 1:2, function(p) sum(diag(var(matrix(p, nrow = n)[, -1, drop = FALSE]))))
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

  # Observed outcomes
  if (object$likelihood %in% c(valid_lhood$binary, valid_lhood$count, "poisson")) {
    y <- object$network$ipd$.r
  } else if (object$likelihood == "ordered") {
    abort("Multinomial models not yet supported.")
    y <- c(t(object$network$ipd$.r[, -1]))
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

