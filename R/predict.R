#' Predictions of absolute effects from NMA models
#'
#' Obtain predictions of absolute effects from NMA models fitted with [nma()].
#' For example, if a model is fitted to binary data with a logit link, predicted
#' outcome probabilities or log odds can be produced.
#'
#' @param object A `stan_nma` object created by [nma()].
#' @param ... Additional arguments, passed to [uniroot()] for regression models
#'   if `baseline_level = "aggregate"`.
#' @param baseline An optional [distr()] distribution for the baseline response
#'   (i.e. intercept), about which to produce absolute effects. If `NULL`,
#'   predictions are produced using the baseline response for each study in the
#'   network with IPD or arm-based AgD.
#'
#'   For regression models, this may be a list of [distr()] distributions of the
#'   same length as the number of studies in `newdata` (possibly named by the
#'   study names, or otherwise in order of appearance in `newdata`).
#'
#'   Use the `baseline_type` and `baseline_level` arguments to specify whether
#'   this distribution is on the response or linear predictor scale, and (for
#'   ML-NMR or models including IPD) whether this applies to an individual at
#'   the reference level of the covariates or over the entire `newdata`
#'   population, respectively. For example, in a model with a logit link with
#'   `baseline_type = "link"`, this would be a distribution for the baseline log
#'   odds of an event. For survival models, `baseline` always corresponds to the
#'   intercept parameters in the linear predictor (i.e. `baseline_type` is
#'   always `"link"`, and `baseline_level` is `"individual"` for IPD NMA or
#'   ML-NMR, and `"aggregate"` for AgD NMA).
#'
#'   Use the `trt_ref` argument to specify which treatment this distribution
#'   applies to.
#' @param newdata Only required if a regression model is fitted and `baseline`
#'   is specified. A data frame of covariate details, for which to produce
#'   predictions. Column names must match variables in the regression model.
#'
#'   If `level = "aggregate"` this should either be a data frame with integration
#'   points as produced by [add_integration()] (one row per study), or a data
#'   frame with individual covariate values (one row per individual) which are
#'   summarised over.
#'
#'   If `level = "individual"` this should be a data frame of individual
#'   covariate values, one row per individual.
#'
#'   If `NULL`, predictions are produced for all studies with IPD and/or
#'   arm-based AgD in the network, depending on the value of `level`.
#' @param study Column of `newdata` which specifies study names or IDs. When not
#'   specified: if `newdata` contains integration points produced by
#'   [add_integration()], studies will be labelled sequentially by row;
#'   otherwise data will be assumed to come from a single study.
#' @param trt_ref Treatment to which the `baseline` response distribution
#'   refers, if `baseline` is specified. By default, the baseline response
#'   distribution will refer to the network reference treatment. Coerced to
#'   character string.
#' @param type Whether to produce predictions on the `"link"` scale (the
#'   default, e.g. log odds) or `"response"` scale (e.g. probabilities).
#'
#'   For survival models, the options are `"survival"` for survival
#'   probabilities (the default), `"hazard"` for hazards, `"cumhaz"` for
#'   cumulative hazards, `"mean"` for mean survival times, `"quantile"` for
#'   quantiles of the survival time distribution, `"median"` for median survival
#'   times (equivalent to `type = "quantile"` with `quantiles = 0.5`), `"rmst"`
#'   for restricted mean survival times, or `"link"` for the linear predictor.
#'   For `type = "survival"`, `"hazard"` or `"cumhaz"`, predictions are given at
#'   the times specified by `times` or at the event/censoring times in the
#'   network if `times = NULL`. For `type = "rmst"`, the restricted time horizon
#'   is specified by `times`, or if `times = NULL` the earliest last follow-up
#'   time amongst the studies in the network is used.
#' @param level The level at which predictions are produced, either
#'   `"aggregate"` (the default), or `"individual"`. If `baseline` is not
#'   specified, predictions are produced for all IPD studies in the network if
#'   `level` is `"individual"` or `"aggregate"`, and for all arm-based AgD
#'   studies in the network if `level` is `"aggregate"`.
#' @param baseline_type When a `baseline` distribution is given, specifies
#'   whether this corresponds to the `"link"` scale (the default, e.g. log odds)
#'   or `"response"` scale (e.g. probabilities). For survival models, `baseline`
#'   always corresponds to the intercept parameters in the linear predictor
#'   (i.e. `baseline_type` is always `"link"`).
#' @param baseline_level When a `baseline` distribution is given, specifies
#'   whether this corresponds to an individual at the reference level of the
#'   covariates (`"individual"`, the default), or from an (unadjusted) average
#'   outcome on the reference treatment in the `newdata` population
#'   (`"aggregate"`). Ignored for AgD NMA, since the only option is
#'   `"aggregate"` in this instance. For survival models, `baseline` always
#'   corresponds to the intercept parameters in the linear predictor (i.e.
#'   `baseline_level` is `"individual"` for IPD NMA or ML-NMR, and `"aggregate"`
#'   for AgD NMA).
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param predictive_distribution Logical, when a random effects model has been
#'   fitted, should the predictive distribution for absolute effects in a new
#'   study be returned? Default `FALSE`.
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @seealso [plot.nma_summary()] for plotting the predictions.
#'
#' @examples ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Predicted log odds of success in each study in the network
#' predict(smk_fit_RE)
#'
#' # Predicted probabilities of success in each study in the network
#' predict(smk_fit_RE, type = "response")
#'
#' # Predicted probabilities in a population with 67 observed events out of 566
#' # individuals on No Intervention, corresponding to a Beta(67, 566 - 67)
#' # distribution on the baseline probability of response, using
#' # `baseline_type = "response"`
#' (smk_pred_RE <- predict(smk_fit_RE,
#'                         baseline = distr(qbeta, 67, 566 - 67),
#'                         baseline_type = "response",
#'                         type = "response"))
#' plot(smk_pred_RE, ref_line = c(0, 1))
#'
#' # Predicted probabilities in a population with a baseline log odds of
#' # response on No Intervention given a Normal distribution with mean -2
#' # and SD 0.13, using `baseline_type = "link"` (the default)
#' # Note: this is approximately equivalent to the above Beta distribution on
#' # the baseline probability
#' (smk_pred_RE2 <- predict(smk_fit_RE,
#'                          baseline = distr(qnorm, mean = -2, sd = 0.13),
#'                          type = "response"))
#' plot(smk_pred_RE2, ref_line = c(0, 1))
#' }
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_mlnmr_example
#' @examples \donttest{
#' # Predicted probabilities of response in each study in the network
#' (pso_pred <- predict(pso_fit, type = "response"))
#' plot(pso_pred, ref_line = c(0, 1))
#'
#' # Predicted probabilites of response in a new target population, with means
#' # and SDs or proportions given by
#' new_agd_int <- data.frame(
#'   bsa_mean = 0.6,
#'   bsa_sd = 0.3,
#'   prevsys = 0.1,
#'   psa = 0.2,
#'   weight_mean = 10,
#'   weight_sd = 1,
#'   durnpso_mean = 3,
#'   durnpso_sd = 1
#' )
#'
#' # We need to add integration points to this data frame of new data
#' # We use the weighted mean correlation matrix computed from the IPD studies
#' new_agd_int <- add_integration(new_agd_int,
#'                                durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
#'                                prevsys = distr(qbern, prob = prevsys),
#'                                bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
#'                                weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
#'                                psa = distr(qbern, prob = psa),
#'                                cor = pso_net$int_cor,
#'                                n_int = 1000)
#'
#' # Predicted probabilities of achieving PASI 75 in this target population, given
#' # a Normal(-1.75, 0.08^2) distribution on the baseline probit-probability of
#' # response on Placebo (at the reference levels of the covariates), are given by
#' (pso_pred_new <- predict(pso_fit,
#'                          type = "response",
#'                          newdata = new_agd_int,
#'                          baseline = distr(qnorm, -1.75, 0.08)))
#' plot(pso_pred_new, ref_line = c(0, 1))
#' }
predict.stan_nma <- function(object, ...,
                             baseline = NULL, newdata = NULL, study = NULL, trt_ref = NULL,
                             type = c("link", "response"),
                             level = c("aggregate", "individual"),
                             baseline_type = c("link", "response"),
                             baseline_level = c("individual", "aggregate"),
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                             predictive_distribution = FALSE,
                             summary = TRUE) {
  # Checks
  if (!inherits(object, "stan_nma")) abort("Expecting a `stan_nma` object, as returned by nma().")

  is_surv <- inherits(object, "stan_nma_surv")  # Survival flag

  if (!is_surv) type <- rlang::arg_match(type)
  level <- rlang::arg_match(level)

  baseline_type <- rlang::arg_match(baseline_type)
  baseline_level <- rlang::arg_match(baseline_level)


  # Get additional survival arguments passed from NextMethod()
  if (is_surv) {
    dlist <- list(...)
    times <- dlist$times
    aux <- dlist$aux
    quantiles <- dlist$quantiles
  }


  # Get network reference treatment
  nrt <- levels(object$network$treatments)[1]

  if (!is.null(trt_ref)) {
    if (is.null(baseline)) {
      # warn("Ignoring `trt_ref` since `baseline` is not given.")
      trt_ref <- nrt
    } else {
      if (length(trt_ref) > 1) abort("`trt_ref` must be length 1.")
      trt_ref <- as.character(trt_ref)
      lvls_trt <- levels(object$network$treatments)
      if (! trt_ref %in% lvls_trt)
        abort(sprintf("`trt_ref` does not match a treatment in the network.\nSuitable values are: %s",
                      ifelse(length(lvls_trt) <= 5,
                             paste0(lvls_trt, collapse = ", "),
                             paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))
    }
  } else {
    # Set trt_ref to network reference treatment if unset
    trt_ref <- nrt
  }

  # Define auxiliary parameters for survival models
  if (is_surv) {
    aux_pars <- switch(object$likelihood,
                       exponential = NULL,
                       `exponential-aft` = NULL,
                       lognormal = "sdlog",
                       gengamma = c("sigma", "k"),
                       mspline = "scoef",
                       pexp = "scoef",
                       "shape")
  } else {
    aux_pars <- NULL
  }

  if (!is_surv || is.null(aux_pars)) {
    if (xor(is.null(newdata), is.null(baseline)) && !is.null(object$regression))
      abort("Specify both `newdata` and `baseline`, or neither.")
  } else {
    if (xor(is.null(newdata), is.null(baseline)) && xor(is.null(newdata), is.null(aux)) && !is.null(object$regression))
      abort("Specify all of `newdata`, `baseline` and `aux`, or none.")
  }

  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) abort("`newdata` is not a data frame.")

    .study <- pull_non_null(newdata, enquo(study))
    if (is.null(.study)) {
      if (inherits(object, "integration_tbl"))
        newdata$.study <- nfactor(paste("New", seq_len(nrow(newdata))))
      else
        newdata$.study <- nfactor("New 1")
    } else {
      check_study(.study)
      newdata <- dplyr::mutate(newdata, .study = nfactor(.study))
    }
  }

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  check_probs(probs)

  # Cannot produce predictions for inconsistency models
  if (object$consistency != "consistency")
    abort(glue::glue("Cannot produce predictions under inconsistency '{object$consistency}' model."))

  # Get NMA formula
  nma_formula <- make_nma_formula(object$regression,
                                  consistency = object$consistency,
                                  classes = !is.null(object$network$classes),
                                  class_interactions = object$class_interactions)


  # Without regression model ---------------------------------------------------
  if (is.null(object$regression)) {

    times <- rlang::eval_tidy(times)

    if (!is.null(baseline)) {
      if (!inherits(baseline, "distr"))
        abort("Baseline response `baseline` should be specified using distr(), or NULL.")
    }

    if (level == "individual")
      abort("Cannot produce individual predictions without a regression model.")

    ## Without baseline specified ----------------------------------------------
    if (is.null(baseline)) {

      if (!has_ipd(object$network) && !has_agd_arm(object$network)) {
        abort("No arm-based data (IPD or AgD) in network. Specify `baseline` to produce predictions of absolute effects.")
      } else {

        # Make design matrix of all studies with baselines, and all treatments
        studies <- forcats::fct_unique(forcats::fct_drop(forcats::fct_c(
          if (has_ipd(object$network)) object$network$ipd$.study else factor(),
          if (has_agd_arm(object$network)) object$network$agd_arm$.study else factor()
          )))
        preddat <- tidyr::expand_grid(.study = studies, .trt = object$network$treatments)

        # Add in .trtclass if defined in network
        if (!is.null(object$network$classes)) {
          preddat$.trtclass <- object$network$classes[as.numeric(preddat$.trt)]
        }

        # Design matrix, just treating all data as AgD arm
        X_list <- make_nma_model_matrix(nma_formula,
                                        dat_agd_arm = preddat,
                                        xbar = object$xbar,
                                        consistency = object$consistency,
                                        classes = !is.null(object$network$classes),
                                        newdata = TRUE)
        X_all <- X_list$X_agd_arm
        rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

        # Get posterior samples
        post <- as.array(object, pars = c("mu", "d"))
        if (predictive_distribution) {
          # For predictive distribution, use delta_new instead of d
          delta_new <- get_delta_new(object)
          post[ , , dimnames(delta_new)[[3]]] <- delta_new
        }

        # Get prediction array
        pred_array <- tcrossprod_mcmc_array(post, X_all)

        # Get aux parameters for survival models
        if (is_surv) {
          if (!is.null(aux)) abort("Specify both `aux` and `baseline` or neither.")

          if (is.null(aux_pars)) {
            aux_array <- NULL
          } else {
            aux_array <- as.array(object, pars = aux_pars)
          }
        }

        # Deal with times argument
        if (type %in% c("survival", "hazard", "cumhaz")) {
          if (is.null(times)) {
            # Take times from network
            surv_all <- dplyr::bind_rows(object$network$ipd,
                                         tidyr::unnest(object$network$agd_arm, cols = ".Surv"))
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Add times vector to preddat
            preddat <- dplyr::left_join(preddat,
                                        dplyr::group_by(surv_all, .data$.study) %>%
                                          dplyr::summarise(.time = list(.data$time)),
                                        by = ".study")

          } else {
            # Use provided vector of times
            if (!is.numeric(times))
              abort("`times` must be a numeric vector of times to predict at.")

            preddat <- dplyr::mutate(preddat, .time = list(times))
          }
        } else if (type == "rmst") {
          if (is.null(times)) {
            # Take time horizon from network
            surv_all <- dplyr::bind_rows(object$network$ipd,
                                         tidyr::unnest(object$network$agd_arm, cols = ".Surv"))
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Time horizon is earliest last follow-up time amongst the studies
            last_times <- surv_all %>%
              dplyr::group_by(.data$.study) %>%
              dplyr::summarise(time = max(time))

            times <- min(last_times$time)
            preddat$.time <- times

          } else {
            # Use provided time
            if (!is.numeric(times) && length(times) > 1)
              abort("`times` must be a scalar numeric value giving the restricted time horizon.")

            preddat$.time <- times
          }
        }

      }
    ## With baseline specified -------------------------------------------------
    } else {

      # Make design matrix of SINGLE study, and all treatments
      preddat <- tibble::tibble(.study = factor("..dummy.."), .trt = object$network$treatments)

      # Add in .trtclass if defined in network
      if (!is.null(object$network$classes)) {
        preddat$.trtclass <- object$network$classes[as.numeric(preddat$.trt)]
      }

      # Design matrix, just treating all data as AgD arm
      X_list <- make_nma_model_matrix(nma_formula,
                                      dat_agd_arm = preddat,
                                      xbar = object$xbar,
                                      consistency = object$consistency,
                                      classes = !is.null(object$network$classes),
                                      newdata = TRUE)
      X_all <- X_list$X_agd_arm
      rownames(X_all) <- paste0("pred[", preddat$.trt, "]")

      # Get posterior samples
      d <- as.array(object, pars = "d")

      # Generate baseline samples
      dim_d <- dim(d)
      dim_mu <- c(dim_d[1:2], 1)
      u <- runif(prod(dim_mu))
      mu <- array(rlang::eval_tidy(rlang::call2(baseline$qfun, p = u, !!! baseline$args)),
                  dim = dim_mu)

      # Convert to linear predictor scale if baseline_type = "response"
      if (baseline_type == "response") {
        mu <- link_fun(mu, link = object$link)
      }

      # Convert to samples on network ref trt if trt_ref given
      if (trt_ref != nrt) {
        mu <- mu - d[ , , paste0("d[", trt_ref, "]"), drop = FALSE]
      }

      # Combine mu and d
      dim_post <- c(dim_d[1:2], dim_d[3] + 1)
      post <- array(NA_real_, dim = dim_post)
      post[ , , 1] <- mu
      if (!predictive_distribution) {
        post[ , , 2:dim_post[3]] <- d
      } else {
        # For predictive distribution, use delta_new instead of d
        post[ , , 2:dim_post[3]] <- get_delta_new(object)
      }

      # Get prediction array
      pred_array <- tcrossprod_mcmc_array(post, X_all)

      # Get aux parameters for survival models
      if (is_surv) {
        if (is.null(aux)) abort("Specify both `aux` and `baseline` or neither.")


        if (is.null(aux_pars)) {
          aux_array <- NULL
        } else {

          if (object$likelihood %in% c("mspline", "pexp")) {
            n_aux <- length(object$basis[[1]])
            aux_names <- paste0(aux_pars, "[..dummy.., ", 1:n_aux, "]")
          } else {
            n_aux <- length(aux_pars)
            aux_names <- paste0(aux_pars, "[..dummy..]")
          }

          dim_aux <- c(dim_mu[1:2], n_aux)
          u <- array(runif(prod(dim_aux)), dim = dim_aux)
          aux_array <- array(NA_real_,
                             dim = dim_aux,
                             dimnames = list(iterations = NULL,
                                             chains = NULL,
                                             parameters = aux_names))

          if (object$likelihood %in% c("mspline", "pexp")) {
            # Generate spline coefficients as a vector
            for (i in 1:dim_aux[1]) {
              for (j in 1:dim_aux[2]) {
                aux_array[i, j, ] <- rlang::eval_tidy(rlang::call2(aux$qfun, p = u[i, j, , drop = TRUE], !!! aux$args))
              }
            }
          } else {
            # All other aux pars generate one by one
            for (i in 1:n_aux) {
              aux_array[, , i] <- rlang::eval_tidy(rlang::call2(aux$qfun, p = u[ , , i, drop = TRUE], !!! aux$args))
            }
          }
        }

        # Check times argument
        if (is.null(times) && type %in% c("survival", "hazard", "cumhaz", "rmst"))
          abort("`times` must be specified when `baseline` and `aux` are provided")

        if (type %in% c("survival", "hazard", "cumhaz") && !is.numeric(times))
          abort("`times` must be a numeric vector of times to predict at.")

        if (type == "rmst" && (!is.numeric(times) && length(times) > 1))
          abort("`times` must be a scalar numeric value giving the restricted time horizon.")

        # Add times vector in to preddat
        if (type %in% c("survival", "hazard", "cumhaz")) {
          preddat <- dplyr::mutate(preddat, .time = list(times))
        } else if (type == "rmst") {
          preddat$.time <- times
        }
      }

    }

    # Get predictions for each category for ordered models
    if (object$likelihood == "ordered") {
      cc <- as.array(object, pars = "cc")
      n_cc <- dim(cc)[3]
      l_cc <- stringr::str_replace(dimnames(cc)[[3]], "^cc\\[(.+)\\]$", "\\1")

      # Apply cutoffs
      d_p <- d_pt <- dim(pred_array)
      d_pt[3] <- d_p[3]*n_cc
      dn_p <- dn_pt <- dimnames(pred_array)
      dn_pt[[3]] <- rep(dn_p[[3]], each = n_cc)
      pred_temp <- array(dim = d_pt, dimnames = dn_pt)

      for (i in 1:d_p[3]) {
        pred_temp[ , , (i-1)*n_cc + 1:n_cc] <- sweep(cc, 1:2, pred_array[ , , i, drop = FALSE],
                                                     FUN = function(x, y) {y - x})
        dn_pt[[3]][(i-1)*n_cc + 1:n_cc] <- paste0(stringr::str_sub(dn_p[[3]][i], start = 1, end = -2),
                                                  ", ", l_cc, "]")
      }
      dimnames(pred_temp) <- dn_pt
      pred_array <- pred_temp
    }

    # Get predictions for survival models
    if (is_surv) {

      pred_temp <- pred_array

      d_p <- dim(pred_temp)
      dn_p <- dimnames(pred_temp)

      if (type %in% c("survival", "hazard", "cumhaz")) {
        d_p[3] <- sum(lengths(preddat$.time))
        dn_p[[3]] <- paste0(rep(stringr::str_sub(dn_p[[3]], start = 1, end = -2), times = lengths(preddat$.time)),
                            ", ", unlist(purrr::map(preddat$.time, seq_along)), "]")
      } else if (type == "quantile") {
        dn_p[[3]] <- paste0(rep(stringr::str_sub(dn_p[[3]], start = 1, end = -2), each = length(quantiles)),
                            ", ", rep(quantiles, times = d_p[3]), "]")
        d_p[3] <- d_p[3]*length(quantiles)
      }

      pred_array <- array(NA_real_, dim = d_p, dimnames = dn_p)

      aux_names <- dimnames(aux_array)[[3]]

      j <- 0
      for (i in 1:nrow(preddat)) {
        if (type %in% c("survival", "hazard", "cumhaz", "rmst")) {
          tt <- preddat$.time[[i]]
          jinc <- length(tt)
        } else if (type == "quantile") {
          tt <- NA
          jinc <- length(quantiles)
        } else {
          tt <- NA
          jinc <- 1
        }
        s <- preddat$.study[i]

        aux_s <- grepl(paste0("[", s, if (object$likelihood %in% c("mspline", "pexp")) "," else "]"),
                       aux_names, fixed = TRUE)

        if (object$likelihood %in% c("mspline", "pexp")) {
          basis <- object$basis[[s]]
        } else {
          basis <- NULL
        }

        pred_array[ , , (j+1):(j+jinc)] <-
          make_surv_predict(eta = pred_temp[ , , i, drop = FALSE],
                            aux = aux_array[ , , aux_s, drop = FALSE],
                            times = tt,
                            quantiles = quantiles,
                            likelihood = object$likelihood,
                            type = type,
                            basis = basis)

        j <- j + jinc
      }
    }

    # Transform predictions if type = "response"
    if (type == "response") {
      pred_array <- inverse_link(pred_array, link = object$link)
    }

    # Produce nma_summary
    if (summary) {
      pred_summary <- summary_mcmc_array(pred_array, probs)

      if (object$likelihood == "ordered") {
        pred_summary <- tibble::add_column(pred_summary,
                                           .trt = rep(preddat$.trt, each = n_cc),
                                           .category = rep(l_cc, times = nrow(preddat)),
                                           .before = 1)

      } else if (object$likelihood %in% valid_lhood$survival) {
        if (type %in% c("survival", "hazard", "cumhaz")) {
          preddat <- tidyr::unnest(preddat, cols = ".time")
          pred_summary <- tibble::add_column(pred_summary,
                                             .trt = preddat$.trt,
                                             .time = preddat$.time,
                                             .before = 1)
        } else if (type == "rmst") {
          pred_summary <- tibble::add_column(pred_summary,
                                             .trt = preddat$.trt,
                                             .time = preddat$.time,
                                             .before = 1)
        } else if (type == "quantile") {
          pred_summary <- tibble::add_column(pred_summary,
                                             .trt = rep(preddat$.trt, each = length(quantiles)),
                                             .quantile = rep(quantiles, times = nrow(preddat)),
                                             .before = 1)
        } else {
          pred_summary <- tibble::add_column(pred_summary,
                                             .trt = preddat$.trt,
                                             .before = 1)
        }

      } else {
        pred_summary <- tibble::add_column(pred_summary,
                                           .trt = preddat$.trt,
                                           .before = 1)
      }

      if (is.null(baseline)) {
        if (object$likelihood == "ordered") {
          pred_summary <- tibble::add_column(pred_summary,
                                             .study = rep(preddat$.study, each = n_cc),
                                             .before = 1)
        } else if (type == "quantile") {
          pred_summary <- tibble::add_column(pred_summary,
                                             .study = rep(preddat$.study, each = length(quantiles)),
                                             .before = 1)
        } else {
          pred_summary <- tibble::add_column(pred_summary,
                                             .study = preddat$.study,
                                             .before = 1)
        }
      }

      out <- list(summary = pred_summary, sims = pred_array)
    } else {
      out <- list(sims = pred_array)
    }

  # With regression model ------------------------------------------------------
  } else {

    if (!is.null(baseline)) {
      if (!(inherits(baseline, "distr") || (rlang::is_list(baseline) && all(purrr::map_lgl(baseline, inherits, what = "distr")))))
        abort("Baseline response `baseline` should be a single distr() specification, a list of distr() specifications, or NULL.")
    }

    ## Without baseline and newdata specified ----------------------------------
    if (is.null(baseline) && is.null(newdata)) {

      times <- rlang::eval_tidy(times)

      # Get data for prediction
      if (level == "individual") {
        if (!has_ipd(object$network))
          abort(paste("No IPD in network to produce individual predictions for.",
                      "  - Specify IPD in `newdata` for which to produce predictions, or",
                      '  - Produce aggregate predictions with level = "aggregate"',
                      sep = "\n"))

        preddat <- get_model_data_columns(object$network$ipd, regression = object$regression)

        if (is_surv) {
          # Deal with times argument
          if (type %in% c("survival", "hazard", "cumhaz")) {
            if (is.null(times)) {
              # Take times from network
              preddat$.time <- get_Surv_data(object$network$ipd$.Surv)$time
            } else {
              abort('Cannot specify `times` with `level = "individual"` and `newdata = NULL`')
            }
          } else if (type == "rmst") {
            if (is.null(times)) {
              # Take time horizon from network
              surv_all <- object$network$ipd %>%
                dplyr::mutate(!!! get_Surv_data(object$network$ipd$.Surv))

              # Time horizon is earliest last follow-up time amongst the studies
              last_times <- surv_all %>%
                dplyr::group_by(.data$.study) %>%
                dplyr::summarise(time = max(time))

              times <- min(last_times$time)
              preddat$.time <- times

            } else {
              # Use provided time
              if (!is.numeric(times) && length(times) > 1)
                abort("`times` must be a scalar numeric value giving the restricted time horizon.")

              preddat$.time <- times
            }
          }

          preddat <- dplyr::group_by(preddat, .data$.study) %>%
            dplyr::mutate(.obs_id = 1:dplyr::n()) %>%
            dplyr::ungroup()
        }

      } else {

        if (!has_ipd(object$network) && !has_agd_arm(object$network)) {
          abort("No arm-based data (IPD or AgD) in network. Specify `baseline` and `newdata` to produce predictions of absolute effects.")
        }

        if ((has_agd_arm(object$network) || has_agd_contrast(object$network)) && !has_agd_sample_size(object$network))
          abort(
            paste("AgD study sample sizes not specified in network, cannot calculate aggregate predictions.",
                  "  - Specify `sample_size` in set_agd_*(), or",
                  "  - Specify covariate values using the `newdata` argument",
                  sep = "\n"))

        # Deal with times argument for type = "rmst"
        if (is_surv && type == "rmst") {
          if (is.null(times)) {
            # Take time horizon from network
            surv_all <- dplyr::bind_rows(object$network$ipd,
                                         tidyr::unnest(object$network$agd_arm, cols = ".Surv"))
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Time horizon is earliest last follow-up time amongst the studies
            last_times <- surv_all %>%
              dplyr::group_by(.data$.study) %>%
              dplyr::summarise(time = max(time))

            times <- min(last_times$time)
          } else {
            # Use provided time
            if (!is.numeric(times) && length(times) > 1)
              abort("`times` must be a scalar numeric value giving the restricted time horizon.")
          }
        } else if (is_surv && type %in% c("survival", "hazard", "cumhaz")) {
          # Checks for survival/hazard/cumhaz times
          if (!is.null(times) && !is.numeric(times))
            abort("`times` must be a numeric vector of times to predict at.")
        }

        if (has_agd_arm(object$network)) {
          if (!is_surv) {
            if (inherits(object$network, "mlnmr_data")) {
              dat_agd_arm <- .unnest_integration(object$network$agd_arm) %>%
                dplyr::mutate(.sample_size = .data$.sample_size / object$network$n_int)
            } else {
              dat_agd_arm <- object$network$agd_arm
            }
          } else {
            if (type %in% c("survival", "hazard", "cumhaz") && is.null(times)) {
                # Use times from network, unnest
                dat_agd_arm <- object$network$agd_arm %>%
                  # Drop duplicated names in outer dataset from .data_orig before unnesting
                  dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(object$network$agd_arm)))),
                                # Reset sample size for weighted mean later
                                .sample_size = 1) %>%
                  tidyr::unnest(cols = c(".Surv", ".data_orig")) %>%
                  dplyr::mutate(!!! get_Surv_data(.$.Surv),
                                .time = .data$time) %>%
                  # Add in ID variable for each observation
                  dplyr::group_by(.data$.study) %>%
                  dplyr::mutate(.obs_id = 1:dplyr::n()) %>%
                  dplyr::ungroup()
            } else {
              # Use provided times
              dat_agd_arm <- object$network$agd_arm %>%
                # Drop duplicated names in outer dataset from .data_orig before unnesting
                # Take only one row of .data_orig (all duplicated anyway)
                dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(object$network$agd_arm)))[1,]),
                              # Use provided time vector
                              .time = if (type %in% c("survival", "hazard", "cumhaz", "rmst")) list(times) else NA,
                              .obs_id = if (type %in% c("survival", "hazard", "cumhaz", "rmst")) list(1:length(times)) else NA,
                              # Reset sample size for weighted mean later
                              .sample_size = 1) %>%
                tidyr::unnest(cols = c(".time", ".obs_id", ".data_orig"))
            }

            # Unnest integration points if present
            if (inherits(object, "stan_mlnmr")) {
              dat_agd_arm <- .unnest_integration(dat_agd_arm) %>%
                dplyr::mutate(.sample_size = .data$.sample_size / object$network$n_int)
            }

            # Drop .Surv column, not needed
            dat_agd_arm <- dplyr::select(dat_agd_arm, -.data$.Surv)
          }

          # Only take necessary columns
          dat_agd_arm <- get_model_data_columns(dat_agd_arm, regression = object$regression, label = "AgD (arm-based)")
        } else {
          dat_agd_arm <- tibble::tibble()
        }

        if (has_ipd(object$network)) {
          dat_ipd <- object$network$ipd

          if (is_surv) {
            # For aggregate survival predictions we average the entire survival
            # curve over the population (i.e. over every individual), so we need
            # to expand out every time for every individual
            if (type %in% c("survival", "hazard", "cumhaz") && is.null(times)) {
              # Use times from network
              ipd_times <- dplyr::mutate(dat_ipd, !!! get_Surv_data(dat_ipd$.Surv), .time = .data$time) %>%
                dplyr::select(".study", ".trt", ".time") %>%
                # Add in ID variable for each observation, needed later for aggregating
                dplyr::group_by(.data$.study) %>%
                dplyr::mutate(.obs_id = 1:dplyr::n()) %>%
                dplyr::group_by(.data$.study, .data$.trt) %>%
                tidyr::nest(.time = .data$.time, .obs_id = .data$.obs_id)

              dat_ipd <- dplyr::left_join(dat_ipd,
                                          ipd_times,
                                          by = c(".study", ".trt")) %>%
                tidyr::unnest(cols = c(".time", ".obs_id"))
            } else {
              # Use provided times
              if (type %in% c("survival", "hazard", "cumhaz", "rmst")) {
                dat_ipd <- dplyr::mutate(dat_ipd,
                                         .time = list(times),
                                         .obs_id = list(1:length(times))) %>%
                  tidyr::unnest(cols = c(".time", ".obs_id"))
              }
            }

            # Drop .Surv column, not needed
            dat_ipd <- dplyr::select(dat_ipd, -.data$.Surv)
          }

          # Only take necessary columns
          dat_ipd <- get_model_data_columns(dat_ipd, regression = object$regression, label = "IPD")

          dat_ipd$.sample_size <- 1
        } else {
          dat_ipd <- tibble::tibble()
        }

        preddat <- dplyr::bind_rows(dat_ipd, dat_agd_arm)
      }

      # Produce predictions on every treatment for each observed arm/individual
      if (packageVersion("dplyr") >= "1.1.1") {
        preddat <- preddat %>%
          dplyr::rename(.trt_old = ".trt") %>%
          dplyr::left_join(tidyr::expand(., .study = .data$.study,
                                            .trt = .data$.trt_old),
                           by = ".study",
                           relationship = "many-to-many")
      } else {
        preddat <- preddat %>%
          dplyr::rename(.trt_old = ".trt") %>%
          dplyr::left_join(tidyr::expand(., .study = .data$.study,
                                         .trt = .data$.trt_old),
                           by = ".study")
      }

      # If producing aggregate-level predictions, output these in factor order
      # Individual-level predictions will be in the order of the input data
      if (level == "aggregate") {
        preddat <- dplyr::arrange(preddat, .data$.study, .data$.trt)
      }

      # Add in .trtclass if defined in network
      if (!is.null(object$network$classes)) {
        preddat$.trtclass <- object$network$classes[as.numeric(preddat$.trt)]
      }

      if (has_agd_contrast(object$network)) {
        dat_agd_contrast <- object$network$agd_contrast
      } else {
        dat_agd_contrast <- tibble::tibble()
      }

      # Design matrix, just treating all data as IPD
      # Contrast data is included just so that the correct columns are excluded
      X_list <- make_nma_model_matrix(nma_formula,
                                      dat_ipd = preddat,
                                      dat_agd_contrast = dat_agd_contrast,
                                      xbar = object$xbar,
                                      consistency = object$consistency,
                                      classes = !is.null(object$network$classes),
                                      newdata = TRUE)
      X_all <- X_list$X_ipd
      rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

      offset_all <- X_list$offset_ipd

      # Get posterior samples
      post <- as.array(object, pars = c("mu", "d", "beta"))

      # Get aux parameters for survival models
      if (is_surv) {
        if (!is.null(aux)) abort("Specify all of `aux`, `baseline`, and `newdata`, or none.")

        if (is.null(aux_pars)) {
          aux_array <- NULL
        } else {
          aux_array <- as.array(object, pars = aux_pars)
        }
      }

    ## With baseline and newdata specified -------------------------------------
    } else {

      if (is_surv) {
        has_int <- inherits(newdata, "integration_tbl")

        newdata <- dplyr::group_by(newdata, .data$.study) %>%
          dplyr::mutate(.obs_id = 1:dplyr::n()) %>%
          dplyr::ungroup()

        # Restore class
        if (has_int) class(newdata) <- c("integration_tbl", class(newdata))
      }

      if (level == "individual") {
        if (!has_ipd(object$network))
          warn("Producing individual predictions from an aggregate-level regression. Interpret with great caution!")

        preddat <- newdata
      } else {
        if (inherits(object, "stan_mlnmr")) {
          if (!inherits(newdata, "integration_tbl")) {
            abort("No integration points found in `newdata`. Specify integration points using add_integration().")
          } else {
            preddat <- .unnest_integration(newdata)
          }
        } else {
          if (has_ipd(object$network) && inherits(newdata, "integration_tbl")) {
            # Allow integration of IPD model over aggregate population
            preddat <- .unnest_integration(newdata)
          } else {
            preddat <- newdata
          }
        }
      }

      # Check all variables are present
      predreg <- get_model_data_columns(preddat, regression = object$regression, label = "`newdata`")

      preddat$.sample_size <- 1

      # Check times argument
      if (is_surv && type %in% c("survival", "hazard", "cumhaz", "rmst")) {
        if (!rlang::is_quosure(times)) times <- rlang::enquo(times)
        if (rlang::quo_is_null(times))
          abort("`times` must be specified when `newdata`, `baseline` and `aux` are provided")


        teval <- rlang::eval_tidy(times, data = preddat)

        if ((rlang::quo_is_symbol(times) || rlang::is_scalar_character(rlang::eval_tidy(times))) && rlang::has_name(preddat, rlang::as_name(times))) {
          # times is a column of newdata
          preddat$.time <- teval

          if (type %in% c("survival", "hazard", "cumhaz") && !is.numeric(teval))
            abort("`times` must be a numeric vector or column of `newdata` giving times to predict at.")

          if (type == "rmst" && !is.numeric(teval))
            abort("`times` must be a scalar numeric value or column of `newdata` giving the restricted time horizon.")

          if (level == "aggregate" && type %in% c("survival", "hazard", "cumhaz")) {
            # Need to average survival curve over all covariate values at every time
            p_times <- dplyr::group_by(preddat, .data$.study) %>%
              dplyr::select(.data$.study, .data$.time, .data$.obs_id) %>%
              tidyr::nest(.time = .data$.time, .obs_id = .data$.obs_id)

            preddat <- dplyr::select(preddat, -.data$.time, -.data$.obs_id) %>%
              dplyr::left_join(p_times, by = ".study") %>%
              tidyr::unnest(cols = c(".time", ".obs_id"))
          }

        } else {
          # times is a vector in global env
          if (type %in% c("survival", "hazard", "cumhaz") && !is.numeric(teval))
            abort("`times` must be a numeric vector or column of `newdata` giving times to predict at.")

          if (type == "rmst" && (!is.numeric(teval) || length(teval) > 1))
            abort("`times` must be a scalar numeric value or column of `newdata` giving the restricted time horizon.")

          # Add times vector in to preddat
          if (type %in% c("survival", "hazard", "cumhaz")) {
            preddat <- dplyr::mutate(preddat, .time = list(teval), .obs_id = list(1:length(teval))) %>%
              tidyr::unnest(cols = c(".time", ".obs_id"))
          } else if (type == "rmst") {
            preddat$.time <- teval
          }
        }
      }

      # Make design matrix of all studies and all treatments
      if (rlang::has_name(preddat, ".trt")) preddat <- dplyr::select(preddat, -".trt")
      if (packageVersion("dplyr") >= "1.1.1") {
        preddat <- dplyr::left_join(preddat,
                                    tidyr::expand(preddat,
                                                  .study = .data$.study,
                                                  .trt = object$network$treatments),
                                    by = ".study",
                                    relationship = "many-to-many")
      } else {
        preddat <- dplyr::left_join(preddat,
                                    tidyr::expand(preddat,
                                                  .study = .data$.study,
                                                  .trt = object$network$treatments),
                                    by = ".study")
      }

      # Add in .trtclass if defined in network
      if (!is.null(object$network$classes)) {
        preddat$.trtclass <- object$network$classes[as.numeric(preddat$.trt)]
      }

      # Design matrix, just treating all data as IPD
      X_list <- make_nma_model_matrix(nma_formula,
                                      dat_ipd = preddat,
                                      xbar = object$xbar,
                                      consistency = object$consistency,
                                      classes = !is.null(object$network$classes),
                                      newdata = TRUE)
      X_all <- X_list$X_ipd
      rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

      offset_all <- X_list$offset_ipd

      # Get posterior samples
      post_temp <- as.array(object, pars = c("d", "beta"))

      # Check baseline
      studies <- unique(preddat$.study)
      n_studies <- length(studies)

      if (!inherits(baseline, "distr")) {
        if (!length(baseline) %in% c(1, n_studies))
          abort(sprintf("`baseline` must be a single distr() distribution, or a list of length %d (number of `newdata` studies)", n_studies))
        if (length(baseline) == 1) {
          baseline <- baseline[[1]]
        } else {
          if (!rlang::is_named(baseline)) {
            names(baseline) <- studies
          } else {
            bl_names <- names(baseline)
            if (dplyr::n_distinct(bl_names) != n_studies)
              abort("`baseline` list names must be distinct study names from `newdata`")
            if (length(bad_bl_names <- setdiff(bl_names, studies)))
              abort(glue::glue("`baseline` list names must match all study names from `newdata`.\n",
                               "Unmatched list names: ",
                               glue::glue_collapse(glue::double_quote(bad_bl_names), sep = ", ", width = 30),
                               ".\n",
                               "Unmatched `newdata` study names: ",
                               glue::glue_collapse(glue::double_quote(setdiff(studies, bl_names)), sep = ", ", width = 30),
                               ".\n"))
          }
        }
      }

      # Generate baseline samples
      dim_post_temp <- dim(post_temp)
      dim_mu <- c(dim_post_temp[1:2], n_studies)
      dimnames_mu <- c(dimnames(post_temp)[1:2], list(parameters = paste0("mu[", levels(studies), "]")))

      if (inherits(baseline, "distr")) {
        u <- runif(prod(dim_mu))
        mu <- array(rlang::eval_tidy(rlang::call2(baseline$qfun, p = u, !!! baseline$args)),
                    dim = dim_mu, dimnames = dimnames_mu)
      } else {
        u <- array(runif(prod(dim_mu)), dim = dim_mu)
        mu <- array(NA_real_, dim = dim_mu, dimnames = dimnames_mu)
        for (s in 1:n_studies) {
          # NOTE: mu must be in *factor order* for later multiplication with design matrix, not observation order
          ss <- levels(studies)[s]
          mu[ , , s] <- array(rlang::eval_tidy(rlang::call2(baseline[[ss]]$qfun, p = u[ , , s], !!! baseline[[ss]]$args)),
                              dim = c(dim_mu[1:2], 1))
        }
      }

      # Convert baseline samples as necessary

      if (!inherits(object, "stan_mlnmr") && !has_ipd(object$network)) {
        # AgD-only regression, ignore baseline_level = "individual"
        # if (baseline_level == "individual")
        #   inform('Setting baseline_level = "aggregate", model intercepts are aggregate level for AgD meta-regression.')

        # Convert to linear predictor scale if baseline_type = "response"
        if (baseline_type == "response") {
          mu <- link_fun(mu, link = object$link)
        }

        # Convert to samples on network ref trt if trt_ref given
        if (trt_ref != nrt) {
          mu <- sweep(mu, 1:2, post_temp[ , , paste0("d[", trt_ref, "]"), drop = FALSE], FUN = "-")
        }
      } else { # ML-NMR or IPD NMR
        if (baseline_level == "individual") {

          # Convert to linear predictor scale if baseline_type = "response"
          if (baseline_type == "response") {
            mu <- link_fun(mu, link = object$link)
          }

          # Convert to samples on network ref trt if trt_ref given
          if (trt_ref != nrt) {
            mu <- sweep(mu, 1:2, post_temp[ , , paste0("d[", trt_ref, "]"), drop = FALSE], FUN = "-")
          }

        } else { # Aggregate baselines

          # Assume that aggregate baselines are *unadjusted*, ie. are crude poolings over reference arm outcomes
          # In this case, we need to marginalise over the natural outcome scale

          if (baseline_type == "link") {
            mu0 <- inverse_link(mu, link = object$link)
          } else {
            mu0 <- mu
          }

          preddat_trt_ref <- dplyr::filter(preddat, .data$.trt == trt_ref)

          # Get posterior samples of betas and d[trt_ref]
          post_beta <- as.array(object, pars = "beta")
          if (trt_ref == nrt) {
            post_d <- 0
          } else {
            post_d <- as.array(object, pars = paste0("d[", trt_ref, "]"))
          }

          # Get design matrix for regression for trt_ref
          X_trt_ref <- X_all[preddat$.trt == trt_ref, , drop = FALSE]
          X_beta_trt_ref <- X_trt_ref[ , !grepl("^(\\.study|\\.trt|\\.contr)[^:]+$", colnames(X_trt_ref)), drop = FALSE]

          if (!is.null(offset_all)) offset_trt_ref <- offset_all[preddat$.trt == trt_ref]

          range_mu <- range(as.array(object, pars = "mu"))

          # Aggregate response on natural scale, use numerical solver
          # Define function to solve for mu
          mu_solve <- function(mu, mu0, post_beta, post_d, X_beta, offset, link) {
            lp <- mu + X_beta %*% post_beta + post_d + offset
            ginv_lp <- inverse_link(lp, link = link)
            return(mu0 - mean(ginv_lp))
          }

          for (s in 1:n_studies) {
            # NOTE: mu must be in *factor order* for later multiplication with design matrix, not observation order

            # Study select
            ss <- preddat_trt_ref$.study == levels(studies)[s]

            s_X_beta <- X_beta_trt_ref[ss, , drop = FALSE]
            if (!is.null(offset_all)) s_offset <- offset_trt_ref[ss]

            for (i_iter in 1:dim_post_temp[1]) {
              for (i_chain in 1:dim_post_temp[2]) {
                rtsolve <- uniroot(mu_solve, interval = range_mu, extendInt = "yes", ...,
                                   mu0 = mu0[i_iter, i_chain, s, drop = TRUE],
                                   post_beta = post_beta[i_iter, i_chain, , drop = TRUE],
                                   post_d = if (trt_ref == nrt) 0 else post_d[i_iter, i_chain, , drop = TRUE],
                                   X_beta = s_X_beta,
                                   offset = if (!is.null(offset_all)) s_offset else 0,
                                   link = object$link)
                mu[i_iter, i_chain, s] <- rtsolve$root
              }
            }
          }
        }
      }

      # Combine mu, d, and beta
      dim_post <- c(dim_post_temp[1:2], dim_mu[3] + dim_post_temp[3])
      dimnames_post <- c(dimnames(post_temp)[1:2], list(parameters = c(dimnames(mu)[[3]], dimnames(post_temp)[[3]])))
      post <- array(NA_real_, dim = dim_post, dimnames = dimnames_post)
      post[ , , 1:dim_mu[3]] <- mu
      post[ , , dim_mu[3] + 1:dim_post_temp[3]] <- post_temp


      # Get aux parameters for survival models
      if (is_surv) {
        if (is.null(aux)) abort("Specify both `aux` and `baseline` or neither.")

        if (is.null(aux_pars)) {
          aux_array <- NULL
        } else {

          # Check aux spec
          if (!inherits(aux, "distr")) {
            if (!length(aux) %in% c(1, n_studies))
              abort(sprintf("`aux` must be a single distr() distribution, or a list of length %d (number of `newdata` studies)", n_studies))
            if (length(aux) == 1) {
              aux <- rep(aux, times = n_studies)
              names(aux) <- studies
            } else {
              if (!rlang::is_named(aux)) {
                names(aux) <- studies
              } else {
                aux_names <- names(aux)
                if (dplyr::n_distinct(aux_names) != n_studies)
                  abort("`aux` list names must be distinct study names from `newdata`")
                if (length(bad_aux_names <- setdiff(aux_names, studies)))
                  abort(glue::glue("`aux` list names must match all study names from `newdata`.\n",
                                   "Unmatched list names: ",
                                   glue::glue_collapse(glue::double_quote(bad_aux_names), sep = ", ", width = 30),
                                   ".\n",
                                   "Unmatched `newdata` study names: ",
                                   glue::glue_collapse(glue::double_quote(setdiff(studies, aux_names)), sep = ", ", width = 30),
                                   ".\n"))
              }
            }
          } else {
            aux <- rep(list(aux), times = n_studies)
            names(aux) <- studies
          }

          if (object$likelihood %in% c("mspline", "pexp")) {
            n_aux <- length(object$basis[[1]])
            aux_names <- paste0(rep(aux_pars, each = n_aux * n_studies), "[", rep(studies, each = n_aux) , rep(1:n_aux, times = n_studies), "]")
          } else {
            n_aux <- length(aux_pars)
            aux_names <- paste0(rep(aux_pars, each = n_studies), "[", rep(studies, times = n_aux) , "]")
          }

          dim_aux <- c(dim_mu[1:2], n_aux * n_studies)
          u <- array(runif(prod(dim_aux)), dim = dim_aux)
          aux_array <- array(NA_real_,
                             dim = dim_aux,
                             dimnames = list(iterations = NULL,
                                             chains = NULL,
                                             parameters = aux_names))

          if (object$likelihood %in% c("mspline", "pexp")) {
            # Generate spline coefficients as a vector
            for (i in 1:dim_aux[1]) {
              for (j in 1:dim_aux[2]) {
                for (s in 1:n_studies) {
                  aux_array[i, j, ((s-1)*n_aux + 1):(s*n_aux)] <- rlang::eval_tidy(rlang::call2(aux$qfun, p = u[i, j, ((s-1)*n_aux + 1):(s*n_aux), drop = TRUE], !!! aux$args))
                }
              }
            }
          } else {
            # All other aux pars generate one by one
            for (s in 1:n_studies) {
              for (i in 1:n_aux) {
                aux_array[, , (s-1)*n_studies + i] <- rlang::eval_tidy(rlang::call2(aux[[studies[s]]]$qfun, p = u[ , , (s-1)*n_studies + i, drop = TRUE], !!! aux[[studies[s]]]$args))
              }
            }
          }
        }
      }
    }

    # For predictive distribution, use delta_new instead of d
    if (predictive_distribution) {
      delta_new <- get_delta_new(object)
      post[ , , dimnames(delta_new)[[3]]] <- delta_new
    }

    # Get cutoffs for ordered models
    if (object$likelihood == "ordered") {
      cc <- as.array(object, pars = "cc")
      n_cc <- dim(cc)[3]
      l_cc <- stringr::str_replace(dimnames(cc)[[3]], "^cc\\[(.+)\\]$", "\\1")
    }

    # Make prediction arrays
    if (is_surv) {
      # Handle survival models separately - produce predictions study by study

      if (level == "individual") {
        if (type %in% c("survival", "hazard", "cumhaz", "rmst")) {
          outdat <- dplyr::select(preddat, .data$.study, .data$.trt, .data$.obs_id, .data$.time)
        } else {
          outdat <- dplyr::select(preddat, .data$.study, .data$.trt, .data$.obs_id)
        }
      } else {
        if (type %in% c("survival", "hazard", "cumhaz")) {
          outdat <- dplyr::distinct(preddat, .data$.study, .data$.trt, .data$.obs_id, .data$.time)
        } else if (type == "rmst") {
          outdat <- dplyr::distinct(preddat, .data$.study, .data$.trt, .data$.time)
        } else {
          outdat <- dplyr::distinct(preddat, .data$.study, .data$.trt)
        }
      }

      if (type == "quantile") {
        outdat <- dplyr::cross_join(outdat,
                                    data.frame(.quantile = quantiles))
      }

      studies <- unique(outdat$.study)
      n_studies <- length(studies)
      treatments <- unique(outdat$.trt)
      n_trt <- length(treatments)

      dim_pred_array <- dim(post)
      dim_pred_array[3] <- nrow(outdat)
      dimnames_pred_array <- dimnames(post)

      if (level == "individual") {
        if (type == "quantile") {
          dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, ", ", outdat$.obs_id, ", ", outdat$.quantile, "]")
        } else {
          dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, ", ", outdat$.obs_id, "]")
        }
      } else if (type %in% c("survival", "hazard", "cumhaz")) {
        dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, ", ", outdat$.obs_id, "]")
      } else if (type == "quantile") {
        dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, ", ", outdat$.quantile, "]")
      } else {
        dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, "]")
      }

      pred_array <- array(NA_real_,
                          dim = dim_pred_array,
                          dimnames = dimnames_pred_array)

      for (s in 1:n_studies) {
        # Study select
        ss <- preddat$.study == studies[s]

        # Aux select
        aux_s <- grepl(paste0("[", studies[s], if (object$likelihood %in% c("mspline", "pexp")) "," else "]"),
                       dimnames(aux_array)[[3]], fixed = TRUE)

        # Linear predictor array for this study
        eta_pred_array <- tcrossprod_mcmc_array(post, X_all[ss, , drop = FALSE])

        if (!is.null(offset_all))
          eta_pred_array <- sweep(eta_pred_array, 3, offset_all[ss], FUN = "+")

        # Basis for mspline models
        if (object$likelihood %in% c("mspline", "pexp")) {
          basis <- object$basis[[s]]
        } else {
          basis <- NULL
        }

        if (type %in% c("survival", "hazard", "cumhaz")) {
          s_time <- preddat$.time[ss]
        } else if (type == "rmst") {
          s_time <- preddat$.time[ss]
          if (length(unique(s_time)) == 1) {
            # Same restriction time for all obs, just take the first (faster)
            s_time <- s_time[1]
          }
        } else {
          s_time <- NULL
        }

        s_pred_array <-
          make_surv_predict(eta = eta_pred_array,
                            aux = aux_array[ , , aux_s, drop = FALSE],
                            times = s_time,
                            quantiles = quantiles,
                            likelihood = object$likelihood,
                            type = type,
                            basis = basis)

        # Aggregate predictions when level = "aggregate"
        if (level == "aggregate") {
          s_preddat <-
            dplyr::select(preddat, ".study", ".trt", ".sample_size", if (type %in% c("survival", "hazard", "cumhaz")) ".obs_id" else NULL) %>%
            dplyr::filter(ss) %>%
            dplyr::mutate(.study = forcats::fct_inorder(forcats::fct_drop(.data$.study)),
                          .trt = forcats::fct_inorder(.data$.trt))

          if (type %in% c("survival", "hazard", "cumhaz")) {
            s_preddat <- dplyr::group_by(s_preddat, .data$.study, .data$.trt, .data$.obs_id)
          } else if (type == "quantile") {
            s_preddat <- dplyr::cross_join(s_preddat,
                                           data.frame(.quantile = quantiles)) %>%
              dplyr::group_by(.data$.study, .data$.trt, .data$.quantile)
          } else {
            s_preddat <- dplyr::group_by(s_preddat, .data$.study, .data$.trt)
          }

          s_preddat <- dplyr::mutate(s_preddat, .weights = .data$.sample_size / sum(.data$.sample_size))

          X_weighted_mean <- Matrix::Matrix(0, ncol = dim(s_pred_array)[3], nrow = dplyr::n_groups(s_preddat))

          X_weighted_mean[cbind(dplyr::group_indices(s_preddat),
                                1:dim(s_pred_array)[3])] <- s_preddat$.weights

          pred_array[ , , outdat$.study == studies[s]] <- tcrossprod_mcmc_array(s_pred_array, X_weighted_mean)
        } else {
          pred_array[ , , outdat$.study == studies[s]] <- s_pred_array
        }

      }

    } else if (level == "individual") {

      # Get prediction array
      pred_array <- tcrossprod_mcmc_array(post, X_all)

      if (!is.null(offset_all))
        pred_array <- sweep(pred_array, 3, offset_all, FUN = "+")

      # Get predictions for each category for ordered models
      if (object$likelihood == "ordered") {
        # Apply cutoffs
        d_p <- d_pt <- dim(pred_array)
        d_pt[3] <- d_p[3]*n_cc
        dn_p <- dn_pt <- dimnames(pred_array)
        dn_pt[[3]] <- rep(dn_p[[3]], each = n_cc)
        pred_temp <- array(dim = d_pt, dimnames = dn_pt)

        for (i in 1:d_p[3]) {
          pred_temp[ , , (i-1)*n_cc + 1:n_cc] <- sweep(cc, 1:2, pred_array[ , , i, drop = FALSE],
                                                       FUN = function(x, y) {y - x})
          dn_pt[[3]][(i-1)*n_cc + 1:n_cc] <- paste0(stringr::str_sub(dn_p[[3]][i], start = 1, end = -2),
                                                    ", ", l_cc, "]")
        }
        dimnames(pred_temp) <- dn_pt
        pred_array <- pred_temp
      }

      # Transform predictions if type = "response"
      if (type == "response") {
        pred_array <- inverse_link(pred_array, link = object$link)
      }

    } else { # Predictions aggregated over each population

      # Produce aggregated predictions study by study - more memory efficient
      outdat <- dplyr::distinct(preddat, .data$.study, .data$.trt)

      if (object$likelihood == "ordered") {
        outdat <- tibble::tibble(.study = rep(outdat$.study, each = n_cc),
                                 .trt = rep(outdat$.trt, each = n_cc),
                                 .cc = rep(l_cc, times = nrow(outdat)))
      }

      studies <- unique(outdat$.study)
      n_studies <- length(studies)
      treatments <- unique(outdat$.trt)
      n_trt <- length(treatments)

      dim_pred_array <- dim(post)
      dim_pred_array[3] <- nrow(outdat)
      dimnames_pred_array <- dimnames(post)
      if (object$likelihood == "ordered") {
        dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, ", ", outdat$.cc, "]")
      } else {
        dimnames_pred_array[[3]] <- paste0("pred[", outdat$.study, ": ", outdat$.trt, "]")
      }

      pred_array <- array(NA_real_,
                          dim = dim_pred_array,
                          dimnames = dimnames_pred_array)

      ss <- vector(length = nrow(outdat))

      for (s in 1:n_studies) {

        # Study select
        ss <- preddat$.study == studies[s]

        # Get prediction array for this study
        s_pred_array <- tcrossprod_mcmc_array(post, X_all[ss, , drop = FALSE])

        if (!is.null(offset_all))
          s_pred_array <- sweep(s_pred_array, 3, offset_all[ss], FUN = "+")

        # Get predictions for each category for ordered models
        if (object$likelihood == "ordered") {
          # Apply cutoffs
          d_p <- d_pt <- dim(s_pred_array)
          d_pt[3] <- d_p[3]*n_cc
          dn_p <- dn_pt <- dimnames(s_pred_array)
          dn_pt[[3]] <- rep(dn_p[[3]], each = n_cc)
          pred_temp <- array(dim = d_pt, dimnames = dn_pt)

          for (i in 1:d_p[3]) {
            pred_temp[ , , (i-1)*n_cc + 1:n_cc] <- sweep(cc, 1:2, s_pred_array[ , , i, drop = FALSE],
                                                         FUN = function(x, y) {y - x})
            dn_pt[[3]][(i-1)*n_cc + 1:n_cc] <- paste0(stringr::str_sub(dn_p[[3]][i], start = 1, end = -2),
                                                      ", ", l_cc, "]")
          }
          dimnames(pred_temp) <- dn_pt
          s_pred_array <- pred_temp
        }

        # Transform predictions if type = "response"
        if (type == "response") {
          s_pred_array <- inverse_link(s_pred_array, link = object$link)
        }

        # Aggregate predictions since level = "aggregate"
        if (object$likelihood == "ordered") {
          s_preddat <- preddat[ss, c(".study", ".trt", ".sample_size")] %>%
            dplyr::slice(rep(seq_len(nrow(.)), each = n_cc)) %>%
            dplyr::mutate(.study = forcats::fct_inorder(forcats::fct_drop(.data$.study)),
                          .trt = forcats::fct_inorder(.data$.trt),
                          .cc = forcats::fct_inorder(rep_len(l_cc, nrow(.)))) %>%
            dplyr::group_by(.data$.study, .data$.trt, .data$.cc) %>%
            dplyr::mutate(.weights = .data$.sample_size / sum(.data$.sample_size))

          X_weighted_mean <- Matrix::Matrix(0, ncol = dim(s_pred_array)[3], nrow = n_trt * n_cc)

          X_weighted_mean[cbind(dplyr::group_indices(s_preddat),
                                1:dim(s_pred_array)[3])] <- s_preddat$.weights
        } else {
          s_preddat <- preddat[ss, c(".study", ".trt", ".sample_size")] %>%
            dplyr::mutate(.study = forcats::fct_inorder(forcats::fct_drop(.data$.study)),
                          .trt = forcats::fct_inorder(.data$.trt)) %>%
            dplyr::group_by(.data$.study, .data$.trt) %>%
            dplyr::mutate(.weights = .data$.sample_size / sum(.data$.sample_size))

          X_weighted_mean <- Matrix::Matrix(0, ncol = dim(s_pred_array)[3], nrow = n_trt)

          X_weighted_mean[cbind(dplyr::group_indices(s_preddat),
                                1:dim(s_pred_array)[3])] <- s_preddat$.weights
        }

        pred_array[ , , outdat$.study == studies[s]] <- tcrossprod_mcmc_array(s_pred_array, X_weighted_mean)

      }

      preddat <- dplyr::distinct(preddat, .data$.study, .data$.trt)
    }

    # Produce nma_summary
    if (summary) {
      pred_summary <- summary_mcmc_array(pred_array, probs)
      if (object$likelihood == "ordered") {
        pred_summary <- tibble::add_column(pred_summary,
                                           .study = rep(preddat$.study, each = n_cc),
                                           .trt = rep(preddat$.trt, each = n_cc),
                                           .category = rep(l_cc, times = nrow(preddat)),
                                           .before = 1)
      } else if (is_surv) {
        pred_summary <- tibble::add_column(pred_summary,
                                           .study = outdat$.study,
                                           .trt = outdat$.trt,
                                           .before = 1)
        if (type %in% c("survival", "hazard", "cumhaz", "rmst")) {
          pred_summary <- tibble::add_column(pred_summary, .time = outdat$.time, .after = ".trt")
        } else if (type == "quantile") {
          pred_summary <- tibble::add_column(pred_summary, .quantile = outdat$.quantile, .after = ".trt")
        }
      } else {
        pred_summary <- tibble::add_column(pred_summary,
                                           .study = preddat$.study,
                                           .trt = preddat$.trt,
                                           .before = 1)
      }
      out <- list(summary = pred_summary, sims = pred_array)
    } else {
      out <- list(sims = pred_array)
    }

  }

  if (summary) {
    attr(out, "xlab") <- "Treatment"
    attr(out, "ylab") <- get_scale_name(likelihood = object$likelihood,
                                        link = object$link,
                                        measure = "absolute",
                                        type = type)
    if (object$likelihood == "ordered") {
      class(out) <- c("ordered_nma_summary", "nma_summary")
    } else if (type %in% c("survival", "hazard", "cumhaz")) {
      class(out) <- c("surv_nma_summary", "nma_summary")
      attr(out, "xlab") <- "Time"
    } else {
      class(out) <- "nma_summary"
    }
  }
  return(out)
}


#' @param times A numeric vector of times to evaluate predictions at.
#'   Alternatively, if `newdata` is specified, `times` can be the name of a
#'   column in `newdata` which contains the times. If `NULL` (the default) then
#'   predictions are made at the event/censoring times from the studies included
#'   in the network. Only used if `type` is `"survival"`, `"hazard"`, `"cumhaz"`
#'   or `"rmst"`.
#' @param aux An optional [distr()] distribution for the auxiliary parameter(s)
#'   in the baseline hazard (e.g. shapes). If `NULL`, predictions are produced
#'   using the parameter estimates for each study in the network with IPD or
#'   arm-based AgD.
#'
#'   For regression models, this may be a list of [distr()] distributions of the
#'   same length as the number of studies in `newdata` (possibly named by the
#'   study names, or otherwise in order of appearance in `newdata`).
#' @param quantiles A numeric vector of quantiles of the survival time
#'   distribution to produce estimates for when `type = "quantile"`.
#'
#' @export
#' @rdname predict.stan_nma
predict.stan_nma_surv <- function(object, times = NULL,
                                  ...,
                                  baseline = NULL,
                                  aux = NULL,
                                  newdata = NULL, study = NULL, trt_ref = NULL,
                                  type = c("survival", "hazard", "cumhaz", "mean", "median", "quantile", "rmst", "link"),
                                  quantiles = c(0.25, 0.5, 0.75),
                                  level = c("aggregate", "individual"),
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                  predictive_distribution = FALSE,
                                  summary = TRUE) {

  type <- rlang::arg_match(type)
  times <- rlang::enquo(times)

  if (!rlang::is_double(quantiles, finite = TRUE) || any(quantiles < 0) || any(quantiles > 1))
    rlang::abort("`quantiles` must be a numeric vector of quantiles between 0 and 1.")

  # Other checks (including times, aux) in predict.stan_nma()
  # Need to pass stan_nma_surv-specific args directly, otherwise these aren't picked up by NextMethod()
  NextMethod(times = times,
             aux = aux,
             type = type,
             quantiles = quantiles,
             baseline_level = "individual",
             baseline_type = "link")
}

#' Produce survival predictions from arrays of linear predictors and auxiliary parameters
#'
#' Designed to work with a single arm at a time
#'
#' @param eta Array of samples from the linear predictor
#' @param aux Array of samples of auxiliary parameters
#' @param times Vector of evaluation times, or time horizon for type = "rmst"
#' @param likelihood Likelihood function to use
#' @param type Type of prediction to create
#' @param quantiles Quantiles for type = "quantile"
#'
#' @noRd
make_surv_predict <- function(eta, aux, times, likelihood,
                              type = c("survival", "hazard", "cumhaz", "mean", "median", "quantile", "rmst", "link"),
                              quantiles = c(0.25, 0.5, 0.75),
                              basis = NULL) {

  if (type == "link") return(eta)

  d_out <- dim(eta)
  dn_out <- dimnames(eta)

  n_eta <- dim(eta)[3]

  if (type %in% c("survival", "hazard", "cumhaz") && n_eta == 1) {
    dn_out[[3]] <- paste0(rep(stringr::str_sub(dn_out[[3]], start = 1, end = -2), each = length(times)),
                          ", ", rep(seq_along(times), times = d_out[3]), "]")
    d_out[3] <- d_out[3] * length(times)
  }

  if (type == "quantile") {
    dn_out[[3]] <- paste0(rep(stringr::str_sub(dn_out[[3]], start = 1, end = -2), each = length(quantiles)),
                          ", ", rep(quantiles, times = d_out[3]), "]")
    d_out[3] <- d_out[3] * length(quantiles)
  }

  out <- array(NA_real_, dim = d_out, dimnames = dn_out)

  if (type %in% c("survival", "hazard", "cumhaz") && n_eta == 1) { # Multiple times, single linear predictor
    for (i in 1:length(times)) {
      out[, , i] <- do.call(surv_predfun(likelihood, type),
                             args = list(times = times[i],
                                         eta = as.vector(eta),
                                         aux = matrix(aux, ncol = dim(aux)[3]),
                                         quantiles = quantiles,
                                         basis = basis))
    }
  } else if (n_eta > 1) { # Single/multiple times, multiple linear predictors
    if (type == "quantile") {
      iinc <- length(quantiles)
      quantiles <- rep(quantiles, each = length(eta[ , , 1]))
    } else {
      iinc <- 1
    }
    for (i in 1:n_eta) {
      ti <- if (length(times) == 1) times else times[i]
      out[ , , ((i-1)*iinc+1):(i*iinc)] <-
        do.call(surv_predfun(likelihood, type),
                args = list(times = ti,
                            eta = as.vector(eta[ , , i]),
                            aux = matrix(aux, ncol = dim(aux)[3]),
                            quantiles = quantiles,
                            basis = basis))
    }
  } else { # Single time, single linear predictor
    if (type == "quantile") quantiles <- rep(quantiles, each = length(eta))
    out <- array(do.call(surv_predfun(likelihood, type),
                         args = list(times = times,
                                     eta = as.vector(eta),
                                     aux = matrix(aux, ncol = dim(aux)[3]),
                                     quantiles = quantiles,
                                     basis = basis)),
                 dim = d_out, dimnames = dn_out)
  }

  return(out)
}

#' Return prediction functions for survival likelihoods
#' @noRd
surv_predfun <- function(likelihood, type) {
  if (likelihood == "exponential") {
    make_predfun(base = "exp", type = type,
                 rate = exp(eta),
                 .ns = list(survival = "stats", median = "stats", quantile = "stats"))

  } else if (likelihood == "weibull") {
    make_predfun(base = "weibullPH", type = type,
                 shape = aux, scale = exp(eta))

  } else if (likelihood == "gompertz") {
    make_predfun(base = "gompertz", type = type,
                 shape = aux, rate = exp(eta))

  } else if (likelihood == "exponential-aft") {
    make_predfun(base = "exp", type = type,
                 rate = exp(-eta),
                 .ns = list(survival = "stats", median = "stats", quantile = "stats"))

  } else if (likelihood == "weibull-aft") {
    make_predfun(base = "weibullPH", type = type,
                 shape = aux, scale = exp(-aux * eta))

  } else if (likelihood == "lognormal") {
    make_predfun(base = "lnorm", type = type,
                 meanlog = eta, sdlog = aux,
                 .ns = list(survival = "stats", median = "stats", quantile = "stats"))

  } else if (likelihood == "loglogistic") {
    make_predfun(base = "llogis", type = type,
                 shape = aux, scale = exp(eta))

  } else if (likelihood == "gamma") {
    make_predfun(base = "gamma", type = type,
                 rate = exp(-eta), shape = aux,
                 .ns = list(survival = "stats", median = "stats", quantile = "stats"))

  } else if (likelihood == "gengamma") {
    make_predfun(base = "gengamma", type = type,
                 mu = eta,
                 sigma = aux[grepl("^sigma\\[", names(aux))],
                 Q = 1 / sqrt(aux[grepl("^k\\[", names(aux))]))

  } else if (likelihood %in% c("mspline", "pexp")) {
    if (likelihood == "mspline" && type == "mean")
      abort("Calculating mean survival time for M-spline models is not supported; basis functions are ill-conditioned past the boundary knots.")
    make_predfun(base = "mspline", type = type,
                 rate = exp(eta), scoef = aux, basis = basis,
                 .ns = NULL)
  }
}

#' Construct prediction functions programmatically
#' @noRd
make_predfun <- function(base, type, ..., .ns = list()) {
  if (rlang::has_name(.ns, type)) .ns <- .ns$type
  else if (!is.null(.ns)) .ns <- "flexsurv"

  if (any(.ns == "flexsurv")) require_pkg("flexsurv")

  fn <- paste0(switch(type,
                      survival = "p",
                      hazard = "h",
                      cumhaz = "H",
                      mean  = "mean_",
                      median = "q",
                      quantile = "q",
                      rmst = "rmst_"),
               base)

  dots <- rlang::enexprs(...)

  if (type == "survival") function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, q = times, lower.tail = FALSE, !!! dots, .ns = .ns))
  else if (type %in% c("hazard", "cumhaz")) function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, x = times, !!! dots, .ns = .ns))
  else if (type == "mean") function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, !!! dots, .ns = .ns))
  else if (type == "median") function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, p = 0.5, !!! dots, .ns = .ns))
  else if (type == "quantile") function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, p = quantiles, !!! dots, .ns = .ns))
  else if (type == "rmst") function(times, eta, aux, quantiles, basis) rlang::eval_tidy(rlang::call2(fn, t = times, !!! dots, .ns = .ns))
}
