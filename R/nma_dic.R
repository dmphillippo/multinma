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
#' @seealso [print.nma_dic()] for printing details, [plot.nma_dic()] for
#'   producing plots of residual deviance contributions.
#'
#' @examples ## Smoking cessation
#' @template ex_smoking_nma_fe_example
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Compare DIC of FE and RE models
#' (smk_dic_FE <- dic(smk_fit_FE))
#' (smk_dic_RE <- dic(smk_fit_RE))   # substantially better fit
#'
#' # Plot residual deviance contributions under RE model
#' plot(smk_dic_RE)
#'
#' # Check for inconsistency using UME model
#' }
#' @template ex_smoking_nma_re_ume_example
#' @examples \donttest{
#' # Compare DIC
#' smk_dic_RE
#' (smk_dic_RE_UME <- dic(smk_fit_RE_UME))  # no difference in fit
#'
#' # Compare residual deviance contributions
#' plot(smk_dic_RE, smk_dic_RE_UME, show_uncertainty = FALSE)
#' }
dic <- function(x, ...) {
  if (!inherits(x, "stan_nma")) abort("Not a `stan_nma` object.")

  net <- x$network

  resdev <- colMeans(as.matrix(x, pars = "resdev"))

  resdev_array <- as.array(x, pars = "resdev")

  if (has_ipd(net)) {
    n_ipd <- nrow(net$ipd)
    resdev_ipd <- resdev[1:n_ipd]
    fitted_ipd <- colMeans(as.matrix(x, pars = "fitted_ipd"))
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(net)) {
    n_agd_arm <- nrow(net$agd_arm)
    resdev_agd_arm <- resdev[n_ipd + (1:n_agd_arm)]
    fitted_agd_arm <- colMeans(as.matrix(x, pars = "fitted_agd_arm"))
  } else {
    n_agd_arm <- 0
  }

  if (has_agd_contrast(net)) {
    # Number of residual deviances is equal to the number of studies
    nr_agd_contrast <- length(unique(net$agd_contrast$.study))
    # Number of fitted values is equal to the number of contrasts
    nf_agd_contrast <- nrow(dplyr::filter(net$agd_contrast, !is.na(.data$.y)))

    resdev_agd_contrast <- resdev[n_ipd + n_agd_arm + (1:nr_agd_contrast)]
    fitted_agd_contrast <- colMeans(as.matrix(x, pars = "fitted_agd_contrast"))
  } else {
    nr_agd_contrast <- nf_agd_contrast <- 0
  }

  has_df <- FALSE

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
      ipd_arm <-  dplyr::group_indices(net$ipd, .data$.study, .data$.trt)

      # Use posterior median for sigma
      ipd_sigma <- apply(as.matrix(x, pars = "sigma"), 2, median)

      resdevfit_ipd <- resdev_ipd
      for (i in seq_along(ipd_y)) {
        resdevfit_ipd[i] <- (ipd_y[i] - fitted_ipd[i])^2 / ipd_sigma[ipd_arm[i]]^2
      }
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

  } else if (x$likelihood == "ordered") {
    if (has_ipd(net)) {
      ipd_r <- net$ipd$.r
      m_fitted_ipd <- matrix(fitted_ipd, nrow = n_ipd)

      resdevfit_ipd <- vector("double", n_ipd)
      for (i in 1:n_ipd) {
        resdevfit_ipd[i] <- 2 * sum((ipd_r[i,] * log(ipd_r[i,] / m_fitted_ipd[i,]))[!is.na(ipd_r[i,]) & ipd_r[i,] > 0])
      }
      leverage_ipd <- resdev_ipd - resdevfit_ipd

      # Degrees of freedom is 1 - number of categories
      df_ipd <- rowSums(!is.na(ipd_r)) - 1
      has_df <- TRUE
    } else {
      leverage_ipd <- NULL
    }

    if (has_agd_arm(net)) {
      agd_arm_r <- net$agd_arm$.r
      m_fitted_agd_arm <- matrix(fitted_agd_arm, nrow = n_agd_arm)

      resdevfit_agd_arm <- vector("double", n_agd_arm)
      for (i in 1:n_agd_arm) {
        resdevfit_agd_arm[i] <- 2 * sum((agd_arm_r[i,] * log(agd_arm_r[i,] / m_fitted_agd_arm[i,]))[!is.na(agd_arm_r[i,]) & agd_arm_r[i,] > 0])
      }
      leverage_agd_arm <- resdev_agd_arm - resdevfit_agd_arm

      # Degrees of freedom is 1 - number of categories
      df_agd_arm <- rowSums(!is.na(agd_arm_r)) - 1
      has_df <- TRUE
    } else {
      leverage_agd_arm <- NULL
    }

    if (has_agd_contrast(net)) {
      df_agd_contrast <- 1
    }

  } else {
    abort(glue::glue("DIC not supported for likelihood of type '{x$likelihood}'."))
  }

  if (has_agd_contrast(net)) {
    # Get covariance structure
    Sigma <- make_Sigma(net$agd_contrast)

    agd_contrast_resdev_dat <-
      net$agd_contrast %>%
        dplyr::filter(!is.na(.data$.y)) %>%
        dplyr::mutate(.fitted = fitted_agd_contrast,
                      .study_inorder = forcats::fct_inorder(forcats::fct_drop(.data$.study))) %>%
        dplyr::group_by(.data$.study_inorder, .data$.study) %>%
        dplyr::summarise(.y = list(.data$.y),
                         fitted = list(.data$.fitted),
                         n_contrast = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(Sigma = Sigma,
                      resdev = resdev_agd_contrast) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(resdevfit = crossprod(.data$.y - .data$fitted,
                                            solve(.data$Sigma,
                                                  .data$.y - .data$fitted)),
                      leverage = .data$resdev - .data$resdevfit)

    leverage_agd_contrast <- agd_contrast_resdev_dat$leverage
  } else {
    leverage_agd_contrast <- NULL
  }

  # Set pointwise contributions
  pw <- list()
  if (has_ipd(net)) {
    pw$ipd <- tibble::tibble(
      .study = net$ipd$.study,
      .trt = net$ipd$.trt,
      resdev = resdev_ipd,
      leverage = leverage_ipd,
      dic = resdev_ipd + leverage_ipd)

    if (has_df) pw$ipd$df <- df_ipd
  } else {
    pw$ipd <- NULL
  }

  if (has_agd_arm(net)) {
    pw$agd_arm <- tibble::tibble(
      .study = net$agd_arm$.study,
      .trt = net$agd_arm$.trt,
      resdev = resdev_agd_arm,
      leverage = leverage_agd_arm,
      dic = resdev_agd_arm + leverage_agd_arm)

    if (has_df) pw$agd_arm$df <- df_agd_arm
  } else {
    pw$agd_arm <- NULL
  }

  if (has_agd_contrast(net)) {
    pw$agd_contrast <- tibble::tibble(
      .study = agd_contrast_resdev_dat$.study,
      n_contrast = agd_contrast_resdev_dat$n_contrast,
      resdev = resdev_agd_contrast,
      leverage = leverage_agd_contrast,
      dic = resdev_agd_contrast + leverage_agd_contrast)

    if (has_df) pw$agd_contrast$df <- df_agd_contrast
  } else {
    pw$agd_contrast <- NULL
  }

  # Calculate pD, DIC
  totresdev <- sum(resdev)
  pd <- sum(leverage_ipd, leverage_agd_arm, leverage_agd_contrast)
  dic <- totresdev + pd

  # Return nma_dic object
  out <- list(dic = dic, pd = pd, resdev = totresdev, pointwise = pw, resdev_array = resdev_array)
  class(out) <- "nma_dic"
  return(out)
}
