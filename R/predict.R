#' Predictions of absolute effects from NMA models
#'
#' Obtain predictions of absolute effects from NMA models fitted with [nma()].
#' For example, if a model is fitted to binary data with a logit link, predicted
#' outcome probabilities or log odds can be produced. For survival models,
#' predictions can be made for survival probabilities, (cumulative) hazards,
#' (restricted) mean survival times, and quantiles including median survival
#' times.
#' When an IPD NMA or ML-NMR model has been fitted, predictions can be
#' produced either at an individual level or at an aggregate level.
#' Aggregate-level predictions are population-average absolute effects; these
#' are marginalised or standardised over a population. For example, average
#' event probabilities from a logistic regression, or marginal (standardised)
#' survival probabilities from a survival model.
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
#'   time amongst the studies in the network is used. When `level =
#'   "aggregate"`, these all correspond to the standardised survival function
#'   (see details).
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
#' @details
#' # Aggregate-level predictions from IPD NMA and ML-NMR models
#' Population-average absolute effects can be produced from IPD NMA and ML-NMR
#' models with `level = "aggregate"`. Predictions are averaged over the target
#' population (i.e. standardised/marginalised), either by (numerical)
#' integration over the joint covariate distribution (for AgD studies in the
#' network for ML-NMR, or AgD `newdata` with integration points created by
#' [add_integration()]), or by averaging predictions for a sample of individuals
#' (for IPD studies in the network for IPD NMA/ML-NMR, or IPD `newdata`).
#'
#' For example, with a binary outcome, the population-average event
#' probabilities on treatment \eqn{k} in study/population \eqn{j} are
#' \deqn{\bar{p}_{jk} = \int_\mathfrak{X} p_{jk}(\mathbf{x}) f_{jk}(\mathbf{x})
#' d\mathbf{x}}
#' for a joint covariate distribution \eqn{f_{jk}(\mathbf{x})} with
#' support \eqn{\mathfrak{X}} or
#' \deqn{\bar{p}_{jk} = \sum_i p_{jk}(\mathbf{x}_i)}
#' for a sample of individuals with covariates \eqn{\mathbf{x}_i}.
#'
#' Population-average absolute predictions follow similarly for other types of
#' outcomes, however for survival outcomes there are specific considerations.
#'
#' ## Standardised survival predictions
#' Different types of population-average survival predictions, often called
#' standardised survival predictions, follow from the **standardised survival
#' function** created by integrating (or equivalently averaging) the
#' individual-level survival functions at each time \eqn{t}:
#' \deqn{\bar{S}_{jk}(t) = \int_\mathfrak{X} S_{jk}(t | \mathbf{x}) f_{jk}(\mathbf{x})
#' d\mathbf{x}}
#' which is itself produced using `type = "survival"`.
#'
#' The **standardised hazard function** corresponding to this standardised
#' survival function is a weighted average of the individual-level hazard
#' functions
#' \deqn{\bar{h}_{jk}(t) = \frac{\int_\mathfrak{X} S_{jk}(t | \mathbf{x}) h_{jk}(t | \mathbf{x}) f_{jk}(\mathbf{x})
#' d\mathbf{x} }{\bar{S}_{jk}(t)}}
#' weighted by the probability of surviving to time \eqn{t}. This is produced
#' using `type = "hazard"`.
#'
#' The corresponding **standardised cumulative hazard function** is
#' \deqn{\bar{H}_{jk}(t) = -\log(\bar{S}_{jk}(t))}
#' and is produced using `type = "cumhaz"`.
#'
#' **Quantiles and medians** of the standardised survival times are found by
#' solving
#' \deqn{\bar{S}_{jk}(t) = 1-\alpha}
#' for the \eqn{\alpha\%} quantile, using numerical root finding. These are
#' produced using `type = "quantile"` or `"median"`.
#'
#' **(Restricted) means** of the standardised survival times are found by
#' integrating
#' \deqn{\mathrm{RMST}_{jk}(t^*) = \int_0^{t^*} \bar{S}_{jk}(t) dt}
#' up to the restricted time horizon \eqn{t^*}, with \eqn{t^*=\infty} for mean
#' standardised survival time. These are produced using `type = "rmst"` or
#' `"mean"`.
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
#'                                n_int = 64)
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
    times_seq <- dlist$times_seq
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
  if (is.null(object$regression) && (is.null(object$aux_by) || length(setdiff(object$aux_by, c(".study", ".trt"))) == 0)) {

    if (is_surv && !is.null(aux_pars) && xor(is.null(baseline), is.null(aux))) {
        abort("Specify both `baseline` and `aux`, or neither")
    }

    if (is_surv) times <- rlang::eval_tidy(times)

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
        if (is.null(object$aux_by) || !".trt" %in% object$aux_by) {
          studies <- forcats::fct_unique(forcats::fct_drop(forcats::fct_c(
            if (has_ipd(object$network)) object$network$ipd$.study else factor(),
            if (has_agd_arm(object$network)) object$network$agd_arm$.study else factor()
            )))
          preddat <- tidyr::expand_grid(.study = studies, .trt = object$network$treatments)
        } else {
          # If aux_by stratifies by .trt too, then we can only predict for observed treatment arms
          preddat <- dplyr::bind_rows(
            if (has_ipd(object$network)) dplyr::distinct(object$network$ipd, .data$.study, .data$.trt) else dplyr::tibble(),
            if (has_agd_arm(object$network)) dplyr::distinct(object$network$agd_arm, .data$.study, .data$.trt) else dplyr::tibble()
          )
        }

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
                                         if (has_agd_arm(object$network)) tidyr::unnest(object$network$agd_arm, cols = ".Surv") else NULL)
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Calculate times sequence if times_seq specified
            if (!is.null(times_seq)) {
              surv_all <- dplyr::group_by(surv_all, .data$.study) %>%
                dplyr::summarise(time = list(seq(from = 0, to = max(.data$time), length.out = times_seq))) %>%
                tidyr::unnest(cols = "time") %>%
                dplyr::ungroup()
            }

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
                                         if (has_agd_arm(object$network)) tidyr::unnest(object$network$agd_arm, cols = ".Surv") else NULL)
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Time horizon is earliest last follow-up time amongst the studies
            last_times <- surv_all %>%
              dplyr::group_by(.data$.study) %>%
              dplyr::summarise(time = max(.data$time))

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

        if (is.null(aux_pars)) {
          aux_array <- NULL
        } else {

          if (object$likelihood %in% c("mspline", "pexp")) {
            abort(glue::glue('Producing predictions for external populations is not currently supported for "{object$likelihood}" models.'))
          } else {
            n_aux <- length(aux_pars)
            aux_names <- paste0(aux_pars, "[..dummy..]")
          }

          if (n_aux > 1) {
            if (!rlang::is_bare_list(aux, n = n_aux) ||
                !setequal(names(aux), aux_pars) ||
                any(purrr::map_lgl(aux, ~!inherits(., "distr"))))
              abort(glue::glue("`aux` must be a named list of distr() specifications for ",
                               glue::glue_collapse(aux_pars, sep = ", ", last = " and "), "."))
          } else {
            if (!inherits(aux, "distr"))
              abort("`aux` must be specified using distr().")
          }

          dim_aux <- c(dim_mu[1:2], n_aux)
          u <- array(runif(prod(dim_aux)), dim = dim_aux)

          if (n_aux == 1) {
            aux_array <- array(rlang::eval_tidy(rlang::call2(aux$qfun, p = u, !!! aux$args)),
                               dim = dim_aux,
                               dimnames = list(iterations = NULL,
                                               chains = NULL,
                                               parameters = aux_names))
          } else {
            aux_array <- array(NA_real_,
                               dim = dim_aux,
                               dimnames = list(iterations = NULL,
                                               chains = NULL,
                                               parameters = aux_names))
            for (i in 1:n_aux) {
              aux_array[, , i] <- rlang::eval_tidy(rlang::call2(aux[[aux_pars[i]]]$qfun, p = u[ , , i, drop = TRUE], !!! aux[[aux_pars[i]]]$args))
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

        s <- get_aux_labels(preddat[i, ], by = object$aux_by)

        aux_s <- grepl(paste0("[", s, if (object$likelihood %in% c("mspline", "pexp")) "," else "]"),
                       aux_names, fixed = TRUE)

        if (object$likelihood %in% c("mspline", "pexp")) {
          basis <- object$basis[[preddat$.study[i]]]
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

    if (!is_surv || is.null(aux_pars)) {
      if (xor(is.null(newdata), is.null(baseline)))
        abort("Specify both `newdata` and `baseline`, or neither.")
    } else {
      if (xor(is.null(newdata), is.null(baseline)) || xor(is.null(newdata), is.null(aux)))
        abort("Specify all of `newdata`, `baseline`, and `aux`, or none.")
    }

    if (!is.null(baseline)) {
      if (!(inherits(baseline, "distr") || (rlang::is_list(baseline) && all(purrr::map_lgl(baseline, inherits, what = "distr")))))
        abort("Baseline response `baseline` should be a single distr() specification, a list of distr() specifications, or NULL.")
    }

    ## Without baseline and newdata specified ----------------------------------
    if (is.null(baseline) && is.null(newdata)) {

      if (is_surv) times <- rlang::eval_tidy(times)

      # Get data for prediction
      if (level == "individual") {
        if (!has_ipd(object$network))
          abort(paste("No IPD in network to produce individual predictions for.",
                      "  - Specify IPD in `newdata` for which to produce predictions, or",
                      '  - Produce aggregate predictions with level = "aggregate"',
                      sep = "\n"))

        preddat <- get_model_data_columns(object$network$ipd, regression = object$regression, keep = object$aux_by)

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
                dplyr::summarise(time = max(.data$time))

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
                                         if (has_agd_arm(object$network)) tidyr::unnest(object$network$agd_arm, cols = ".Surv") else NULL
                                         )
            surv_all <- dplyr::mutate(surv_all, !!! get_Surv_data(surv_all$.Surv))

            # Time horizon is earliest last follow-up time amongst the studies
            last_times <- surv_all %>%
              dplyr::group_by(.data$.study) %>%
              dplyr::summarise(time = max(.data$time))

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
              if (is.null(times_seq)) {
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
                # Calculate times sequence if times_seq specified
                agd_times <- object$network$agd_arm %>%
                  tidyr::unnest(cols = ".Surv") %>%
                  dplyr::mutate(!!! get_Surv_data(.$.Surv), .time = .data$time) %>%
                  dplyr::group_by(.data$.study) %>%
                  dplyr::summarise(.time = list(seq(from = 0, to = max(.data$.time), length.out = times_seq)),
                                   .obs_id = list(1:times_seq))

                dat_agd_arm <- object$network$agd_arm %>%
                  # Drop duplicated names in outer dataset from .data_orig before unnesting
                  # Take only one row of .data_orig (all duplicated anyway)
                  dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(object$network$agd_arm)))[1,]),
                                # Reset sample size for weighted mean later
                                .sample_size = 1) %>%
                  dplyr::left_join(agd_times, by = ".study") %>%
                  tidyr::unnest(cols = c(".time", ".obs_id", ".data_orig"))
              }
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
            dat_agd_arm <- dplyr::select(dat_agd_arm, -".Surv")
          }

          # Only take necessary columns
          dat_agd_arm <- get_model_data_columns(dat_agd_arm, regression = object$regression, keep = object$aux_by, label = "AgD (arm-based)")
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
              if (is.null(times_seq)) {
                # Use times from network
                ipd_times <- dplyr::mutate(dat_ipd, !!! get_Surv_data(dat_ipd$.Surv), .time = .data$time) %>%
                  dplyr::select(".study", ".trt", ".time") %>%
                  # Add in ID variable for each observation, needed later for aggregating
                  dplyr::group_by(.data$.study) %>%
                  dplyr::mutate(.obs_id = 1:dplyr::n()) %>%
                  dplyr::group_by(.data$.study, .data$.trt) %>%
                  tidyr::nest(.time = ".time", .obs_id = ".obs_id")

                dat_ipd <- dplyr::left_join(dat_ipd,
                                            ipd_times,
                                            by = c(".study", ".trt")) %>%
                  tidyr::unnest(cols = c(".time", ".obs_id"))

              } else {
                # Calculate times sequence if times_seq specified
                ipd_times <- dplyr::mutate(dat_ipd, !!! get_Surv_data(dat_ipd$.Surv), .time = .data$time) %>%
                  dplyr::group_by(.data$.study) %>%
                  dplyr::summarise(.time = list(seq(from = 0, to = max(.data$.time), length.out = times_seq)),
                                   .obs_id = list(1:times_seq)) %>%
                  dplyr::ungroup()

                dat_ipd <- dplyr::left_join(dat_ipd,
                                            ipd_times,
                                            by = ".study") %>%
                  tidyr::unnest(cols = c(".time", ".obs_id"))
              }

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
            dat_ipd <- dplyr::select(dat_ipd, -".Surv")
          }

          # Only take necessary columns
          dat_ipd <- get_model_data_columns(dat_ipd, regression = object$regression, keep = object$aux_by, label = "IPD")

          dat_ipd$.sample_size <- 1
        } else {
          dat_ipd <- tibble::tibble()
        }

        preddat <- dplyr::bind_rows(dat_ipd, dat_agd_arm)
      }

      # Produce predictions on every treatment for each observed arm/individual
      if (is.null(object$aux_by) || ! ".trt" %in% object$aux_by) {
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
            inform(paste0("No integration points found in `newdata`, averaging over individuals provided in `newdata`.\n",
                          "To set up integration points, use add_integration()."))
            preddat <- newdata
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
      predreg <- get_model_data_columns(preddat, regression = object$regression, keep = object$aux_by, label = "`newdata`")

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
              dplyr::select(".study", ".time", ".obs_id") %>%
              tidyr::nest(.time = ".time", .obs_id = ".obs_id")

            preddat <- dplyr::select(preddat, -".time", -".obs_id") %>%
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

        if (is.null(aux_pars)) {
          aux_array <- NULL
        } else {

          # Check aux spec
          n_aux <- length(aux_pars)

          if (n_aux == 1) {
            if (!inherits(aux, "distr") && length(aux) != n_studies)
              abort(sprintf("`aux` must be a single distr() specification, or a list of length %d (number of `newdata` studies)", n_studies))
            if (inherits(aux, "distr")) {
              aux <- rep(list(aux), times = n_studies)
              names(aux) <- studies
            } else {
              if (any(purrr::map_lgl(aux, ~!inherits(., "distr"))))
                  abort(sprintf("`aux` must be a single distr() specification, or a list of length %d (number of `newdata` studies)", n_studies))
              if (!rlang::is_named(aux)) {
                names(aux) <- studies
              } else {
                aux_names <- names(aux)
                if (!setequal(aux_names, studies))
                  abort(glue::glue("`aux` list names must match all study names from `newdata`.\n",
                                   "Unmatched list names: ",
                                   glue::glue_collapse(glue::double_quote(setdiff(aux_names, studies)), sep = ", ", width = 30),
                                   ".\n",
                                   "Unmatched `newdata` study names: ",
                                   glue::glue_collapse(glue::double_quote(setdiff(studies, aux_names)), sep = ", ", width = 30),
                                   ".\n"))
              }
            }
          } else {
            aux_names <- names(aux)
            if (!rlang::is_bare_list(aux) ||
                !length(aux) %in% c(n_aux, n_studies) ||
                (!setequal(aux_names, aux_pars) && !setequal(aux_names, aux_pars)) ||
                any(purrr::map_lgl(purrr::list_flatten(aux), ~!inherits(., "distr")))) {
              abort(glue::glue("`aux` must be a single named list of distr() specifications for {glue::glue_collapse(aux_pars, sep = ', ', last = ' and ')}, ",
                               "or a list of length {n_studies} (number of `newdata` studies) of such lists."))
            }

            if (setequal(aux_names, aux_pars)) {
              aux <- rep(list(aux), times = n_studies)
              names(aux) <- studies
            } else if (!rlang::is_named(aux)) {
              names(aux) <- studies
            }
          }

          if (object$likelihood %in% c("mspline", "pexp")) {
            abort(glue::glue('Producing predictions for external populations is not currently supported for "{object$likelihood}" models.'))
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

          if (n_aux == 1) {
            for (s in 1:n_studies) {
              aux_array[, , s] <- rlang::eval_tidy(rlang::call2(aux[[studies[s]]]$qfun, p = u[ , , s, drop = TRUE], !!! aux[[studies[s]]]$args))
            }
          } else {
            for (s in 1:n_studies) {
              for (i in 1:n_aux) {
                aux_array[, , (s-1)*n_studies + i] <- rlang::eval_tidy(rlang::call2(aux[[studies[s]]][[i]]$qfun, p = u[ , , (s-1)*n_studies + i, drop = TRUE], !!! aux[[studies[s]]][[i]]$args))
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
          outdat <- dplyr::select(preddat, ".study", ".trt", ".obs_id", ".time")
        } else {
          outdat <- dplyr::select(preddat, ".study", ".trt", ".obs_id")
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

        # Basis for mspline models
        if (object$likelihood %in% c("mspline", "pexp")) {
          basis <- object$basis[[s]]
        } else {
          basis <- NULL
        }

        for (trt in 1:n_trt) {

          # Select corresponding rows
          ss <- preddat$.study == studies[s] & preddat$.trt == treatments[trt]

          # Linear predictor array for this study
          eta_pred_array <- tcrossprod_mcmc_array(post, X_all[ss, , drop = FALSE])

          if (!is.null(offset_all))
            eta_pred_array <- sweep(eta_pred_array, 3, offset_all[ss], FUN = "+")

          # Aux select
          aux_l <- get_aux_labels(preddat[ss, ], by = object$aux_by)
          aux_id <- get_aux_id(preddat[ss, ], by = object$aux_by)

          if (length(setdiff(object$aux_by, c(".study", ".trt"))) == 0) {
            aux_s <- grepl(paste0("\\[(", paste(aux_l, collapse = "|"), if (object$likelihood %in% c("mspline", "pexp")) ")," else ")\\]"),
                           dimnames(aux_array)[[3]])
            aux_array_s <- aux_array[ , , aux_s, drop = FALSE]

          } else {
            # Stratified aux pars within arm, need to expand these out over the individuals/integration points

            if (object$likelihood %in% c("mspline", "pexp")) {
              aux_array_s <- array(NA_real_, dim = c(dim(eta_pred_array), length(basis)))
              for (i in 1:length(basis)) {
                aux_s <- grep(paste0("\\[(", paste(aux_l, collapse = "|"), "), ", i, "\\]"),
                              dimnames(aux_array)[[3]])
                aux_array_s[ , , , i] <- aux_array[ , , aux_s[aux_id]]
              }
            } else if (object$likelihood == "gengamma") {
              aux_array_s <- array(NA_real_, dim = c(dim(eta_pred_array), 2),
                                   dimnames = list(iterations = NULL, chains = NULL, parameters = NULL, aux = c("sigma[]", "k[]")))
              sigma_s <- grep(paste0("^sigma\\[(", paste(aux_l, collapse = "|"), ")", "\\]"),
                            dimnames(aux_array)[[3]])
              k_s <- grep(paste0("^k\\[(", paste(aux_l, collapse = "|"), ")", "\\]"),
                              dimnames(aux_array)[[3]])
              aux_array_s[ , , , "sigma[]"] <- aux_array[ , , sigma_s[aux_id]]
              aux_array_s[ , , , "k[]"] <- aux_array[ , , k_s[aux_id]]

            } else {
              aux_s <- grep(paste0("\\[(", paste(aux_l, collapse = "|"), ")\\]"),
                            dimnames(aux_array)[[3]])
              aux_array_s <- aux_array[ , , aux_s[aux_id], drop = FALSE]
            }
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

          # Aggregate predictions when level = "aggregate"
          if (level == "aggregate") {
            s_preddat <-
              dplyr::select(preddat, ".study", ".trt", ".sample_size", if (type %in% c("survival", "hazard", "cumhaz")) ".obs_id" else NULL) %>%
              dplyr::filter(ss) %>%
              dplyr::mutate(.study = forcats::fct_inorder(forcats::fct_drop(.data$.study)),
                            .trt = forcats::fct_inorder(forcats::fct_drop(.data$.trt)))

            if (type %in% c("survival", "hazard", "cumhaz")) {
              s_preddat <- dplyr::group_by(s_preddat, .data$.study, .data$.trt, .data$.obs_id)
            } else {
              s_preddat <- dplyr::group_by(s_preddat, .data$.study, .data$.trt)
            }

            s_preddat <- dplyr::mutate(s_preddat, .weights = .data$.sample_size / sum(.data$.sample_size))

            pred_array[ , , outdat$.study == studies[s] & outdat$.trt == treatments[trt]] <-
              make_agsurv_predict(eta = eta_pred_array,
                                  aux = aux_array_s,
                                  times = s_time,
                                  quantiles = quantiles,
                                  likelihood = object$likelihood,
                                  type = type,
                                  basis = basis,
                                  weights = s_preddat$.weights,
                                  id = dplyr::group_indices(s_preddat))


          } else { # Individual predictions
            s_pred_array <-
              make_surv_predict(eta = eta_pred_array,
                                aux = aux_array_s,
                                times = s_time,
                                quantiles = quantiles,
                                likelihood = object$likelihood,
                                type = type,
                                basis = basis)

            pred_array[ , , outdat$.study == studies[s] & outdat$.trt == treatments[trt]] <- s_pred_array
          }

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
#'   in the network (or according to `times_seq`). Only used if `type` is
#'   `"survival"`, `"hazard"`, `"cumhaz"` or `"rmst"`.
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
#' @param times_seq A positive integer, when specified evaluate predictions at
#'   this many evenly-spaced event times between 0 and the latest follow-up time
#'   in each study, instead of every observed event/censoring time. Only used if
#'   `newdata = NULL` and `type` is `"survival"`, `"hazard"` or `"cumhaz"`. This
#'   can be useful for plotting survival or (cumulative) hazard curves, where
#'   prediction at every observed even/censoring time is unnecessary and can be
#'   slow. When a call from within `plot()` is detected, e.g. like
#'   `plot(predict(fit, type = "survival"))`, `times_seq` will default to 50.
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
                                  times_seq = NULL,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                  predictive_distribution = FALSE,
                                  summary = TRUE) {

  type <- rlang::arg_match(type)
  times <- rlang::enquo(times)

  if (!rlang::is_double(quantiles, finite = TRUE) || any(quantiles < 0) || any(quantiles > 1))
    abort("`quantiles` must be a numeric vector of quantiles between 0 and 1.")

  if (!is.null(times_seq) && (!rlang::is_integerish(times_seq, n = 1, finite = TRUE) || times_seq <= 0))
    abort("`times_seq` must be a single positive integer.")

  # Set times_seq by default if called within plot()
  if (is.null(times_seq) && type %in% c("survival", "hazard", "cumhaz") && deparse(sys.call(-2)[[1]]) == "plot")
    times_seq <- 50

  # Other checks (including times, aux) in predict.stan_nma()
  # Need to pass stan_nma_surv-specific args directly, otherwise these aren't picked up by NextMethod()
  NextMethod(times = times,
             aux = aux,
             type = type,
             quantiles = quantiles,
             times_seq = times_seq,
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
#' @param basis M-spline basis
#'
#' @noRd
make_surv_predict <- function(eta, aux, times, likelihood,
                              type = c("survival", "hazard", "cumhaz", "mean", "median", "quantile", "rmst", "link"),
                              quantiles = c(0.25, 0.5, 0.75),
                              basis = NULL) {

  if (type == "link") return(eta)

  # Check for times beyond mspline boundary knots
  if (likelihood == "mspline") {
    lower <- attr(basis, "Boundary.knots")[1]
    upper <- attr(basis, "Boundary.knots")[2]
    if (! type %in% c("median", "quantile") && (type == "mean" || any(times < lower) || any(times > upper))) warn("Evaluating M-spline at times beyond the boundary knots.")
  }

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
    if (!is.null(aux)) aux <- matrix(aux, ncol = dim(aux)[3], dimnames = list(NULL, dimnames(aux)[[3]]))
    for (i in 1:length(times)) {
      out[, , i] <- do.call(surv_predfuns[[likelihood]][[type]],
                            args = list(times = times[i],
                                        eta = as.vector(eta),
                                        aux = aux,
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
      auxi <- NULL
      if (!is.null(aux)) {
        if (length(dim(aux)) == 4) {
          auxi <- matrix(aux[ , , i, ], ncol = dim(aux)[4], dimnames = list(NULL, dimnames(aux)[[4]]))
        } else if (likelihood %in% c("gengamma", "mspline", "pexp")) {
          auxi <- matrix(aux, ncol = dim(aux)[3], dimnames = list(NULL, dimnames(aux)[[3]]))
        } else if (dim(aux)[3] == 1) {
          auxi <- as.vector(aux)
        } else {
          auxi <- as.vector(aux[ , , i])
        }
      }
      out[ , , ((i-1)*iinc+1):(i*iinc)] <-
        do.call(surv_predfuns[[likelihood]][[type]],
                args = list(times = ti,
                            eta = as.vector(eta[ , , i]),
                            aux = auxi,
                            quantiles = quantiles,
                            basis = basis))
    }
  } else { # Single time, single linear predictor
    if (!is.null(aux)) aux <- matrix(aux, ncol = dim(aux)[3], dimnames = list(NULL, dimnames(aux)[[3]]))
    if (type == "quantile") quantiles <- rep(quantiles, each = length(eta))
    out <- array(do.call(surv_predfuns[[likelihood]][[type]],
                         args = list(times = times,
                                     eta = as.vector(eta),
                                     aux = aux,
                                     quantiles = quantiles,
                                     basis = basis)),
                 dim = d_out, dimnames = dn_out)
  }

  # Check for times beyond mspline boundary knots
  if (likelihood == "mspline" && type %in% c("median", "quantile")) {
    if (any(out < lower) || any(out > upper)) warn("Evaluating M-spline at times beyond the boundary knots.")
  }

  return(out)
}

#' Produce aggregate (standardised) survival predictions from arrays of linear predictors and auxiliary parameters
#'
#' Designed to work with a single arm at a time
#'
#' @param eta Array of samples from the linear predictor
#' @param aux Array of samples of auxiliary parameters
#' @param times Vector of evaluation times, or time horizon for type = "rmst"
#' @param likelihood Likelihood function to use
#' @param type Type of prediction to create
#' @param quantiles Quantiles for type = "quantile"
#' @param basis M-spline basis
#' @param weights Weights for weighted mean over the population
#' @param id Time ID for survival/hazard predictions
#'
#' @noRd
make_agsurv_predict <- function(eta, aux, times, likelihood,
                                type = c("survival", "hazard", "cumhaz", "mean", "median", "quantile", "rmst", "link"),
                                quantiles = c(0.25, 0.5, 0.75),
                                basis = NULL, weights, id) {

  if (type %in% c("survival", "cumhaz", "link")) {

    X_weighted_mean <- Matrix::Matrix(0, ncol = dim(eta)[3], nrow = max(id))

    X_weighted_mean[cbind(id, 1:dim(eta)[3])] <- weights

    pred_array <-
      make_surv_predict(eta = eta,
                        aux = aux,
                        times = times,
                        quantiles = quantiles,
                        likelihood = likelihood,
                        type = if (type == "cumhaz") "survival" else type,
                        basis = basis)

    out <- tcrossprod_mcmc_array(pred_array, X_weighted_mean)

    if (type == "cumhaz") out <- -log(out)

  } else if (type == "hazard") {

    X_weighted_mean <- Matrix::Matrix(0, ncol = dim(eta)[3], nrow = max(id))

    X_weighted_mean[cbind(id, 1:dim(eta)[3])] <- weights

    S_array <-
      make_surv_predict(eta = eta,
                        aux = aux,
                        times = times,
                        quantiles = quantiles,
                        likelihood = likelihood,
                        type = "survival",
                        basis = basis)

    h_array <-
      make_surv_predict(eta = eta,
                        aux = aux,
                        times = times,
                        quantiles = quantiles,
                        likelihood = likelihood,
                        type = "hazard",
                        basis = basis)

    Shbar <- tcrossprod_mcmc_array(S_array * h_array, X_weighted_mean)
    Sbar <- tcrossprod_mcmc_array(S_array, X_weighted_mean)
    out <- Shbar / Sbar

    # Fix up edge case with all S=0
    S0 <- which(Shbar == 0 & Sbar == 0)
    out[S0] <- 0

  } else if (type %in% c("median", "quantile")) {

    if (type == "median") quantiles <- 0.5

    out <- array(NA_real_, dim = c(dim(eta)[1:2], length(quantiles)))

    for (i in 1:length(quantiles)) {
      out[ , , i] <- quantile_Sbar(p = quantiles[i],
                           eta = eta,
                           weights = weights,
                           aux = aux,
                           likelihood = likelihood,
                           basis = basis)
    }

  } else if (type %in% c("mean", "rmst")) {

    if (type == "mean") times <- Inf

    out <- array(NA_real_, dim = c(dim(eta)[1:2], length(times)))

    for (i in 1:dim(out)[1]) for (j in 1:dim(out)[2]) {
      out[i, j, ] <- rmst_Sbar(times = times,
                               eta = eta[i, j, ],
                               weights = weights,
                               aux = if (length(dim(aux)) == 3) matrix(aux[i, j, ], ncol = dim(aux)[3], dimnames = list(NULL, dimnames(aux)[[3]])) else aux[i, j, , ],
                               likelihood = likelihood,
                               basis = basis)
    }
  }

  # Check for times beyond mspline boundary knots
  if (likelihood == "mspline") {
    lower <- attr(basis, "Boundary.knots")[1]
    upper <- attr(basis, "Boundary.knots")[2]
    if (type == "mean" || any(out < lower) || any(out > upper)) warn("Evaluating M-spline at times beyond the boundary knots.")
  }

  return(out)
}

#' Restricted mean survival times for marginal survival curves
#' Given linear predictor vector evaluated over the population and associated weights
#' @noRd
rmst_Sbar <- function(times, eta, weights, aux, likelihood, basis, start = 0) {
  require_pkg("flexsurv")

  flexsurv::rmst_generic(pSbar, t = times, start = start,
                         eta = eta, weights = weights, aux = aux,
                         likelihood = likelihood, basis = basis,
                         scalarargs = c("basis", "likelihood", "eta", "aux", "weights"))
}

pSbar <- function(times, eta, weights, aux, likelihood, basis) {
  Sbar <- numeric(length(times))
  for (i in 1:length(times)) {
    S <- do.call(surv_predfuns[[likelihood]][["survival"]],
                       args = list(times = times[i],
                                   eta = eta,
                                   aux = aux,
                                   basis = basis))

    Sbar[i] <- weighted.mean(S, weights)
  }
  return(1 - Sbar)
}


#' Quantile function for marginal survival curves
#' Given linear predictor array evaluated over the population and associated weights
#' @noRd
quantile_Sbar <- function(p, eta, weights, aux, likelihood, basis) {
  # require_pkg("flexsurv")
  #
  # q <- flexsurv::qgeneric(Sbar, p = rep(p, each = prod(dim(eta)[1:2])),
  #                         eta = eta, weights = weights, aux = aux,
  #                         likelihood = likelihood, basis = basis,
  #                         scalarargs = c("basis", "likelihood", "eta", "aux", "weights"),
  #                         lbound = 0)
  #
  # return(q)

  # Faster to call rstpm2::vuniroot directly

  require_pkg("rstpm2")

  upper <- if (!is.null(basis)) attr(basis, "Boundary.knots")[2] else 1

  rstpm2::vuniroot(qSbar,
                   lower = 0,
                   upper = upper,
                   extendInt = "downX",
                   n = prod(dim(eta)[1:2]),
                   p = 1 - p,
                   eta = eta,
                   weights = weights,
                   aux = aux,
                   likelihood = likelihood,
                   basis = basis,
                   tol = .Machine$double.eps)$root

}

qSbar <- function(times, p, eta, weights, aux, likelihood, basis) {
  S <- array(NA_real_, dim = dim(eta))
  auxi <- NULL
  for (i in 1:dim(eta)[3]) {
    if (!is.null(aux)) {
      if (length(dim(aux)) == 4) {
        auxi <- matrix(aux[ , , i, ], ncol = dim(aux)[4], dimnames = list(NULL, dimnames(aux)[[4]]))
      } else if (likelihood %in% c("gengamma", "mspline", "pexp")) {
        auxi <- matrix(aux, ncol = dim(aux)[3], dimnames = list(NULL, dimnames(aux)[[3]]))
      } else if (dim(aux)[3] == 1) {
        auxi <- as.vector(aux)
      } else {
        auxi <- as.vector(aux[ , , i])
      }
    }
    S[ , , i] <- do.call(surv_predfuns[[likelihood]][["survival"]],
                         args = list(times = times,
                                     eta = as.vector(eta[ , , i]),
                                     aux = auxi,
                                     basis = basis))
  }

  Sbar <- apply(S, MARGIN = 1:2, FUN = weighted.mean, w = weights)

  Sbar - p
}


#' Construct prediction functions programmatically
#' @noRd
make_predfun <- function(base, type, ..., .ns = list()) {

  dots <- rlang::enexprs(...)

  if (!is.null(.ns)) {
    .ns <- purrr::list_modify(list(survival = "flexsurv", hazard = "flexsurv", cumhaz = "flexsurv",
                                   mean = "flexsurv", quantile = "flexsurv", rmst = "flexsurv"),
                              !!! .ns)
  }

  list(
    survival = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$survival == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("p", base), q = quote(times), lower.tail = FALSE, !!! dots, .ns = .ns$survival)
      })),
    hazard = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$hazard == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("h", base), x = quote(times), !!! dots, .ns = .ns$hazard)
      })),
    cumhaz = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$cumhaz == "flexsurv") quote(require_pkg("flexsurv"))
        !! rlang::call2(paste0("H", base), x = quote(times), !!! dots, .ns = .ns$cumhaz)
      })),
    mean = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$mean == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("mean_", base), !!! dots, .ns = .ns$mean)
      })),
    median = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$quantile == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("q", base), p = 0.5, !!! dots, .ns = .ns$quantile)
      })),
    quantile = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$quantile == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("q", base), p = quote(quantiles), !!! dots, .ns = .ns$quantile)
      })),
    rmst = rlang::new_function(rlang::pairlist2(times=, eta=, aux=, quantiles=, basis=), rlang::expr({
      !! if (!is.null(.ns) && .ns$rmst == "flexsurv") quote(require_pkg("flexsurv"))
      !! rlang::call2(paste0("rmst_", base), t = quote(times), !!! dots, .ns = .ns$rmst)
      }))
  )

}

surv_predfuns <- list(
  exponential = make_predfun(base = "exp",
                             rate = exp(eta),
                             .ns = list(survival = "stats", median = "stats", quantile = "stats")),

  weibull = make_predfun(base = "weibullPH",
                         shape = aux, scale = exp(eta)),

  gompertz = make_predfun(base = "gompertz",
                          shape = aux, rate = exp(eta)),

  `exponential-aft` = make_predfun(base = "exp",
                                   rate = exp(-eta),
                                   .ns = list(survival = "stats", median = "stats", quantile = "stats")),

  `weibull-aft` = make_predfun(base = "weibullPH",
                               shape = aux, scale = exp(-aux * eta)),

  lognormal = make_predfun(base = "lnorm",
                           meanlog = eta, sdlog = aux,
                           .ns = list(survival = "stats", median = "stats", quantile = "stats")),

  loglogistic = make_predfun(base = "llogis",
                             shape = aux, scale = exp(eta)),

  gamma = make_predfun(base = "gamma",
                       rate = exp(-eta), shape = aux,
                       .ns = list(survival = "stats", median = "stats", quantile = "stats")),

  gengamma = make_predfun(base = "gengamma",
                          mu = eta,
                          sigma = aux[, grepl("^sigma\\[", colnames(aux))],
                          Q = 1 / sqrt(aux[, grepl("^k\\[", colnames(aux))])),

  mspline = make_predfun(base = "mspline",
                         rate = exp(eta), scoef = aux, basis = basis,
                         .ns = NULL),

  pexp = make_predfun(base = "mspline",
                      rate = exp(eta), scoef = aux, basis = basis,
                      .ns = NULL)
)
