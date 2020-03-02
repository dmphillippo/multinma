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

  resdev <- colMeans(as.matrix(x, pars = "resdev"))
  fitted <- colMeans(as.matrix(x, pars = "fitted"))

  resdev_array <- as.array(x, pars = "resdev")
  dn_resdev_array <- dimnames(resdev_array)

  if (has_ipd(net)) {
    n_ipd <- nrow(net$ipd)
    resdev_ipd <- resdev[1:n_ipd]
    fitted_ipd <- fitted[1:n_ipd]

    dn_resdev_array[[3]][1:n_ipd] <-
      paste0("resdev[", make_data_labels(net$ipd$.study, net$ipd$.trt), "]")
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(net)) {
    n_agd_arm <- nrow(net$agd_arm)
    resdev_agd_arm <- resdev[n_ipd + (1:n_agd_arm)]
    fitted_agd_arm <- fitted[n_ipd + (1:n_agd_arm)]

    dn_resdev_array[[3]][n_ipd + (1:n_agd_arm)] <-
      paste0("resdev[", make_data_labels(net$agd_arm$.study, net$agd_arm$.trt), "]")

  } else {
    n_agd_arm <- 0
  }

  if (has_agd_contrast(net)) {
    # Number of residual deviances is equal to the number of studies
    nr_agd_contrast <- length(unique(net$agd_contrast$.study))
    # Number of fitted values is equal to the number of contrasts
    nf_agd_contrast <- nrow(dplyr::filter(net$agd_contrast, !is.na(.data$.y)))

    resdev_agd_contrast <- resdev[n_ipd + n_agd_arm + (1:nr_agd_contrast)]
    fitted_agd_contrast <- fitted[n_ipd + n_agd_arm + (1:nf_agd_contrast)]

    # dn_resdev_array is fixed below for agd_contrast (one resdev per study)

  } else {
    nr_agd_contrast <- nf_agd_contrast <- 0
  }

  # Apply formatted dimnames to resdev_array
  dimnames(resdev_array) <- dn_resdev_array


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

  } else {
    abort(glue::glue("DIC not supported for likelihood of type '{x$likelihood}'."))
  }

  if (has_agd_contrast(net)) {
    # Get covariance structure
    Sigma <- make_Sigma(net$agd_contrast)

    agd_contrast_resdev_dat <-
      net$agd_contrast %>%
        dplyr::filter(!is.na(.data$.y)) %>%
        dplyr::mutate(.fitted = fitted_agd_contrast) %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::summarise(.y = list(.data$.y),
                         fitted = list(.data$.fitted),
                         n_contrast = dplyr::n()) %>%
        dplyr::mutate(Sigma = Sigma, resdev = resdev_agd_contrast) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(resdevfit = tcrossprod(crossprod(.data$.y - .data$fitted,
                                                       solve(.data$Sigma)),
                                             .data$.y - .data$fitted),
                      leverage = .data$resdev - .data$resdevfit)

    leverage_agd_contrast <- agd_contrast_resdev_dat$leverage

    dn_resdev_array[[3]][n_ipd + n_agd_arm + (1:nf_agd_contrast)] <-
      paste0(agd_contrast_resdev_dat$.study, " (",
             agd_contrast_resdev_dat$n_contrast + 1, " arms)")

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

#' Make labels for data points
#'
#' @param study Study vector, coerced to character
#' @param trt Treatment vector, coerced to character
#' @param trt_b Baseline treatment vector if contrast-based, coerced to character
#'
#' @return Character vector of labels
#'
#' @noRd
make_data_labels <- function(study, trt, trt_b = NA) {
  tibble::tibble(study = study, trt = trt, trt_b = trt_b) %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(grp_n = dplyr::n(),
                  grp_id = 1:dplyr::n(),
                  vs_trt_b = dplyr::if_else(is.na(trt_b),
                                            NA_character_,
                                            paste0(" vs. ", trt_b)),
                  label = as.character(
                    dplyr::if_else(.data$grp_n > 1,
                      glue::glue("{study}: {trt}{vs_trt_b}, {grp_id}", .na = ""),
                      glue::glue("{study}: {trt}{vs_trt_b}", .na = "")))
    ) %>%
    dplyr::pull(.data$label)
}
