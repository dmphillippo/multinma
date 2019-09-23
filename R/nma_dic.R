#' Deviance Information Criterion (DIC)
#'
#' Calculate the DIC for a model fitted using the [nma()] function.
#'
#' @param x A fitted model object, inheriting class [stan_nma]
#' @param ... Other arguments (not used)
#'
#' @return A [nma_dic] object.
#' @export
#'
#' @examples
dic <- function(x, ...) {
  if (!inherits(x, "stan_nma")) abort("Not a `stan_nma` object.")

  net <- x$network
  sf <- as.stanfit(x)

  resdev <- colMeans(as.matrix(sf, pars = "resdev"))
  fitted <- colMeans(as.matrix(sf, pars = "fitted"))

  if (has_ipd(net)) {
    n_ipd <- nrow(net$ipd)
    resdev_ipd <- resdev[1:n_ipd]
    fitted_ipd <- fitted[1:n_ipd]
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(net)) {
    n_agd_arm <- nrow(net$agd_arm)
    resdev_agd_arm <- resdev[n_ipd + (1:n_agd_arm)]
    fitted_agd_arm <- fitted[n_ipd + (1:n_agd_arm)]
  } else {
    n_agd_arm <- 0
  }

  if (has_agd_contrast(net)) {
    n_agd_contrast <- nrow(net$agd_contrast)
    resdev_agd_contrast <- resdev[n_ipd + n_agd_arm + (1:n_agd_contrast)]
    fitted_agd_contrast <- fitted[n_ipd + n_agd_arm + (1:n_agd_contrast)]
  } else {
    n_agd_contrast <- 0
  }


  if (x$likelihood %in% c("bernoulli", "bernoulli2", "binomial", "binomial2")) {
    if (has_ipd(net)) {
      ipd_r <- net$ipd$.r
      resdevfit_ipd <- 2 * ifelse(ipd_r == 1,
                                  ipd_r * log(ipd_r / fitted_ipd),
                                  (1 - ipd_r) * log((1 - ipd_r) / (1 - fitted_ipd)))
      leverage_ipd <- resdev_ipd - resdevfit_ipd
    } else {
      leverage_ipd <- NULL
    }

    if (has_agd_arm(net)) {
      agd_arm_r <- net$agd_arm$.r
      agd_arm_n <- net$agd_arm$.n
      resdevfit_agd_arm <- 2 * (ifelse(agd_arm_r > 0,
                                       agd_arm_r * log(agd_arm_r / fitted_agd_arm),
                                       0) +
                                ifelse(agd_arm_n - agd_arm_r > 0,
                                       (agd_arm_n - agd_arm_r) *
                                         log((agd_arm_n - agd_arm_r) / (agd_arm_n - fitted_agd_arm)),
                                       0))
      leverage_agd_arm <- resdev_agd_arm - resdevfit_agd_arm
    } else {
      leverage_agd_arm <- NULL
    }

  } else if (x$likelihood == "normal") {
    if (has_ipd(net)) {
      ipd_y <- net$ipd$.y
      # Use posterior median for sigma
      ipd_sigma <- median(as.matrix(sf, pars = "sigma"))
      resdevfit_ipd <- (ipd_y - fitted_ipd)^2 / ipd_sigma^2
      leverage_ipd <- resdev_ipd - resdevfit_ipd
    } else {
      leverage_ipd <- NULL
    }

    if (has_agd_arm(net)) {
      agd_arm_y <- net$agd_arm$.y
      agd_arm_se <- net$agd_arm$.se
      resdevfit_agd_arm <- (agd_arm_y - fitted_agd_arm)^2 / agd_arm_se^2
      leverage_agd_arm <- resdev_agd_arm - resdevfit_agd_arm
    } else {
      leverage_agd_arm <- NULL
    }

  } else if (x$likelihood == "poisson") {
    if (has_ipd(net)) {
      ipd_r <- net$ipd$.r
      resdevfit_ipd <- 2 * ((fitted_ipd - ipd_r) +
                              ifelse(ipd_r > 0, ipd_r * log(ipd_r / fitted_ipd), 0))
      leverage_ipd <- resdev_ipd - resdevfit_ipd
    } else {
      leverage_ipd <- NULL
    }

    if (has_agd_arm(net)) {
      agd_arm_r <- net$agd_arm$.r
      resdevfit_agd_arm <- 2 * ((fitted_agd_arm - agd_arm_r) +
                              ifelse(agd_arm_r > 0, agd_arm_r * log(agd_arm_r / fitted_agd_arm), 0))
      leverage_agd_arm <- resdev_agd_arm - resdevfit_agd_arm
    } else {
      leverage_agd_arm <- NULL
    }

  } else {
    abort(glue::glue("DIC not supported for likelihood of type '{x$likelihood}'."))
  }

  if (has_agd_contrast(net)) {
    agd_contrast_y <- net$agd_contrast$.y
    agd_contrast_se <- net$agd_contrast$.se
    resdevfit_agd_contrast <- (agd_contrast_y - fitted_agd_contrast)^2 / agd_contrast_se^2
    leverage_agd_contrast <- resdev_agd_contrast - resdevfit_agd_contrast
  } else {
    leverage_agd_contrast <- NULL
  }

  # Set pointwise contributions
  pw <- list()
  if (has_ipd(net)) {
    pw$ipd <- tibble::tibble(resdev = resdev_ipd,
                             leverage = leverage_ipd,
                             dic = resdev_ipd + leverage_ipd)
  } else {
    pw$ipd <- NULL
  }

  if (has_agd_arm(net)) {
    pw$agd_arm <- tibble::tibble(resdev = resdev_agd_arm,
                             leverage = leverage_agd_arm,
                             dic = resdev_agd_arm + leverage_agd_arm)
  } else {
    pw$agd_arm <- NULL
  }

  if (has_agd_contrast(net)) {
    pw$agd_contrast <- tibble::tibble(resdev = resdev_agd_contrast,
                                 leverage = leverage_agd_contrast,
                                 dic = resdev_agd_contrast + leverage_agd_contrast)
  } else {
    pw$agd_contrast <- NULL
  }

  # Calculate pD, DIC
  dbar <- sum(resdev)
  pd <- sum(leverage_ipd, leverage_agd_arm, leverage_agd_contrast)
  dic <- dbar + pd

  # Return nma_dic object
  out <- list(dic = dic, pd = pd, dbar = dbar, pointwise = pw)
  class(out) <- "nma_dic"
  return(out)
}
