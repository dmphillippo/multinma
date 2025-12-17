#' Deviance Information Criterion (DIC)
#'
#' Calculate the DIC for a model fitted using the [nma()] function.
#'
#' @param x A fitted model object, inheriting class [stan_nma]
#' @param penalty The method for estimating the effective number of parameters,
#'   used to penalise model fit in the DIC. Either `"pD"` (the default), or
#'   `"pV"`. For survival likelihoods only `"pV"` is currently available.
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
dic <- function(x, penalty = c("pD", "pV"), ...) {
  if (!inherits(x, "stan_nma")) abort("Not a `stan_nma` object.")

  penalty <- rlang::arg_match(penalty)

  if (x$likelihood %in% valid_lhood$survival) penalty <- "pV"

  net <- x$network

  resdev <- colMeans(as.matrix(x, pars = "resdev"))

  resdev_array <- as.array(x, pars = "resdev")

  if (has_ipd(net)) {
    n_ipd <- nrow(net$ipd)
    resdev_ipd <- resdev[1:n_ipd]
    fitted_ipd <- if (!x$likelihood %in% valid_lhood$survival) colMeans(as.matrix(x, pars = "fitted_ipd")) else NULL
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(net)) {
    if (!x$likelihood %in% valid_lhood$survival) {
      n_agd_arm <- nrow(net$agd_arm)
      resdev_agd_arm <- resdev[n_ipd + (1:n_agd_arm)]
      fitted_agd_arm <- colMeans(as.matrix(x, pars = "fitted_agd_arm"))
    } else {
      n_agd_arm <- sum(net$agd_arm$.sample_size)
      resdev_agd_arm <- resdev[n_ipd + (1:n_agd_arm)]
      fitted_agd_arm <- NULL
    }
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
  df_ipd <- df_agd_arm <- df_agd_contrast <- 1

  if (x$likelihood %in% c("bernoulli", "bernoulli2", "binomial", "binomial2") && penalty != "pV") {
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

  } else if (x$likelihood == "normal" && penalty != "pV") {
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

  } else if (x$likelihood == "poisson" && penalty != "pV") {
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
      fitted_ipd <- matrix(fitted_ipd, nrow = n_ipd)

      if (penalty != "pV") {
        resdevfit_ipd <- vector("double", n_ipd)
        for (i in 1:n_ipd) {
          resdevfit_ipd[i] <- 2 * sum((ipd_r[i,] * log(ipd_r[i,] / fitted_ipd[i,]))[!is.na(ipd_r[i,]) & ipd_r[i,] > 0])
        }
        leverage_ipd <- resdev_ipd - resdevfit_ipd
      }

      # Degrees of freedom is 1 - number of categories
      df_ipd <- rowSums(!is.na(ipd_r)) - 1
      has_df <- TRUE
    } else {
      leverage_ipd <- NULL
    }

    if (has_agd_arm(net)) {
      agd_arm_r <- net$agd_arm$.r
      fitted_agd_arm <- matrix(fitted_agd_arm, nrow = n_agd_arm)

      if (penalty != "pV") {
        resdevfit_agd_arm <- vector("double", n_agd_arm)
        for (i in 1:n_agd_arm) {
          resdevfit_agd_arm[i] <- 2 * sum((agd_arm_r[i,] * log(agd_arm_r[i,] / fitted_agd_arm[i,]))[!is.na(agd_arm_r[i,]) & agd_arm_r[i,] > 0])
        }
        leverage_agd_arm <- resdev_agd_arm - resdevfit_agd_arm
      }

      # Degrees of freedom is 1 - number of categories
      df_agd_arm <- rowSums(!is.na(agd_arm_r)) - 1
      has_df <- TRUE
    } else {
      leverage_agd_arm <- NULL
    }

  } else if (x$likelihood %in% valid_lhood$survival) {
    # Nothing to do, pV only

  } else if (penalty != "pV") {
    abort(glue::glue("DIC not supported for likelihood of type '{x$likelihood}'."))
  }

  if (has_agd_contrast(net)) {
    # Get covariance structure
    Sigma <- make_Sigma(net$agd_contrast)

    agd_contrast_resdev_dat <-
      net$agd_contrast %>%
        dplyr::filter(!is.na(.data$.y)) %>%
        dplyr::mutate(.fitted = fitted_agd_contrast,
                      .observed = get_outcome_variables(., net$outcome$agd_contrast)[[1]],
                      .study_inorder = forcats::fct_inorder(forcats::fct_drop(.data$.study))) %>%
        dplyr::group_by(.data$.study_inorder, .data$.study) %>%
        dplyr::summarise(.y = list(.data$.y),
                         fitted = list(.data$.fitted),
                         observed = list(.data$.observed),
                         n_contrast = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(Sigma = Sigma,
                      resdev = resdev_agd_contrast) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(resdevfit = drop(crossprod(.data$.y - .data$fitted,
                                            solve(.data$Sigma,
                                                  .data$.y - .data$fitted))),
                      leverage = .data$resdev - .data$resdevfit)

    leverage_agd_contrast <- agd_contrast_resdev_dat$leverage

    # df for multi-arm trials
    if (any(agd_contrast_resdev_dat$n_contrast > 1)) {
      has_df <- TRUE
      df_agd_contrast <- agd_contrast_resdev_dat$n_contrast
    }
  } else {
    leverage_agd_contrast <- NULL
  }

  # Get pointwise contributions for pV
  if (penalty == "pV") {
    leverage <- colSums(var(matrix(resdev_array, ncol = dim(resdev_array)[3])))/2
    if (has_ipd(net)) leverage_ipd <- leverage[1:n_ipd]
    if (has_agd_arm(net)) leverage_agd_arm <- leverage[n_ipd + (1:n_agd_arm)]
    if (has_agd_contrast(net)) leverage_agd_contrast <- leverage[n_ipd + n_agd_arm + (1:nr_agd_contrast)]
  }

  # Set pointwise contributions
  pw <- list()
  if (has_ipd(net)) {
    pw$ipd <- tibble::tibble(
      .study = net$ipd$.study,
      .trt = net$ipd$.trt,
      resdev = resdev_ipd,
      leverage = leverage_ipd,
      dic = resdev_ipd + leverage_ipd,
      fitted = fitted_ipd,
      observed = get_outcome_variables(net$ipd, net$outcome$ipd)[[1]])

    if (has_df) pw$ipd$df <- df_ipd
  } else {
    pw$ipd <- NULL
  }

  if (has_agd_arm(net)) {
    aa <- net$agd_arm
    if (x$likelihood %in% valid_lhood$survival) aa <- tidyr::unnest(aa, cols = ".Surv")

    pw$agd_arm <- tibble::tibble(
      .study = aa$.study,
      .trt = aa$.trt,
      resdev = resdev_agd_arm,
      leverage = leverage_agd_arm,
      dic = resdev_agd_arm + leverage_agd_arm,
      fitted = fitted_agd_arm,
      observed = get_outcome_variables(aa, net$outcome$agd_arm)[[1]])

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
      dic = resdev_agd_contrast + leverage_agd_contrast,
      fitted = agd_contrast_resdev_dat$fitted,
      observed = agd_contrast_resdev_dat$observed)

    if (has_df) pw$agd_contrast$df <- df_agd_contrast
  } else {
    pw$agd_contrast <- NULL
  }

  # Calculate pD, DIC
  totresdev <- sum(resdev)

  if (penalty == "pD") {
    pd <- sum(leverage_ipd, leverage_agd_arm, leverage_agd_contrast)
    dic <- totresdev + pd
    pv <- NULL
  } else if (penalty == "pV") {
    pv <- var(rowSums(as.matrix(x, pars = "resdev"))) / 2
    dic <- totresdev + pv
    pd <- NULL
  }

  # Return nma_dic object
  out <- list(dic = dic, pd = pd, pv = pv, resdev = totresdev, pointwise = pw, resdev_array = resdev_array)
  class(out) <- "nma_dic"
  attr(out, "penalty") <- penalty
  return(out)
}
