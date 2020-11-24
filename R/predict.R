#' Predictions of absolute effects from NMA models
#'
#' Obtain predictions of absolute effects from NMA models fitted with [nma()].
#' For example, if a model is fitted to binary data with a logit link, predicted
#' outcome probabilities or log odds can be produced.
#'
#' @param object A `stan_nma` object created by [nma()].
#' @param ... Additional arguments (not used).
#' @param baseline An optional [distr()] distribution for the baseline response
#'   (i.e. intercept) on the linear predictor scale, about which to produce
#'   absolute effects. For example, in a model with a logit link, this would be
#'   a distribution for the baseline log odds of an event. If `NULL`,
#'   predictions are produced using the baseline response for each study in the
#'   network with IPD or arm-based AgD.
#' @param newdata Only required if a regression model is fitted and `baseline`
#'   is specified. A data frame of covariate details, for which to produce
#'   predictions. Column names must match variables in the regression model.
#'
#'   If `type = "aggregate"` this should either be a data frame with integration
#'   points as produced by [add_integration()] (one row per study), or a data
#'   frame with individual covariate values (one row per individual) which are
#'   summarised over.
#'
#'   If `type = "individual"` this should be a data frame of individual
#'   covariate values, one row per individual.
#'
#'   If `NULL`, predictions are produced for all studies with IPD and/or
#'   arm-based AgD in the network, depending on the value of `type`.
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
#' @param level The level at which predictions are produced, either
#'   `"aggregate"` (the default), or `"individual"`. If `baseline` is not
#'   specified, predictions are produced for all IPD studies in the network if
#'   `type` is `"individual"` or `"aggregate"`, and for all arm-based AgD
#'   studies in the network if `type` is `"aggregate"`.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
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
#' # Predicted probabilities in a population with a baseline log odds of
#' # response on No Intervention given a Normal distribution with mean -2
#' # and SD 0.15
#' (smk_pred_RE <- predict(smk_fit_RE,
#'                         baseline = distr(qnorm, mean = -2, sd = 0.15),
#'                         type = "response"))
#' plot(smk_pred_RE, ref_line = c(0, 1))
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
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                             summary = TRUE) {
  # Checks
  if (!inherits(object, "stan_nma")) abort("Expecting a `stan_nma` object, as returned by nma().")

  type <- rlang::arg_match(type)
  level <- rlang::arg_match(level)

  if (!is.null(baseline)) {
    if (!inherits(baseline, "distr"))
      abort("Baseline response `baseline` should be specified using distr(), or NULL.")
  }

  if (!is.null(trt_ref)) {
    if (is.null(baseline)) {
      # warn("Ignoring `trt_ref` since `baseline` is not given.")
      trt_ref <- NULL
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
  }

  if (xor(is.null(newdata), is.null(baseline)) && !is.null(object$regression))
    abort("Specify both `newdata` and `baseline`, or neither.")

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
    abort(glue::glue("Cannot produce predictions under inconsistency '{x$consistency}' model."))

  # Get network reference treatment
  nrt <- levels(object$network$treatments)[1]

  # Get NMA formula
  nma_formula <- make_nma_formula(object$regression,
                                  consistency = object$consistency,
                                  classes = !is.null(object$network$classes),
                                  class_interactions = object$class_interactions)

  # Without regression model
  if (is.null(object$regression)) {

    if (level == "individual")
      abort("Cannot produce individual predictions without a regression model.")

    # Without baseline specified
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
                                        classes = !is.null(object$network$classes))
        X_all <- X_list$X_agd_arm
        rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

        # Get posterior samples
        post <- as.array(object, pars = c("mu", "d"))

        # Get prediction array
        pred_array <- tcrossprod_mcmc_array(post, X_all)

      }
    # With baseline specified
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
                                      classes = !is.null(object$network$classes))
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

      # Convert to samples on network ref trt if trt_ref given
      if (!is.null(trt_ref) && trt_ref != nrt) {
        mu <- mu - d[ , , paste0("d[", trt_ref, "]"), drop = FALSE]
      }

      # Combine mu and d
      dim_post <- c(dim_d[1:2], dim_d[3] + 1)
      post <- array(NA_real_, dim = dim_post)
      post[ , , 1] <- mu
      post[ , , 2:dim_post[3]] <- d

      # Get prediction array
      pred_array <- tcrossprod_mcmc_array(post, X_all)

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

    # Transform predictions if type = "response"
    if (type == "response") {
      pred_array <- inverse_link(pred_array, link = object$link)
    }

    # Produce nma_summary
    if (summary) {
      pred_summary <- summary_mcmc_array(pred_array, probs)
      if (is.null(baseline)) {
        if (object$likelihood == "ordered") {
          pred_summary <- tibble::add_column(pred_summary,
                                             .study = rep(preddat$.study, each = n_cc),
                                             .before = 1)
        } else {
          pred_summary <- tibble::add_column(pred_summary, .study = preddat$.study, .before = 1)
        }
      }
      out <- list(summary = pred_summary, sims = pred_array)
    } else {
      out <- list(sims = pred_array)
    }

  # With regression model
  } else {

    # Without baseline and newdata specified
    if (is.null(baseline) && is.null(newdata)) {

      # Get data for prediction
      if (level == "individual") {
        if (!has_ipd(object$network))
          abort(paste("No IPD in network to produce individual predictions for.",
                      "  - Specify IPD in `newdata` for which to produce predictions, or",
                      '  - Produce aggregate predictions with level = "aggregate"',
                      sep = "\n"))

        preddat <- object$network$ipd

      } else {
        if ((has_agd_arm(object$network) || has_agd_contrast(object$network)) && !has_agd_sample_size(object$network))
          abort(
            paste("AgD study sample sizes not specified in network, cannot calculate aggregate predictions.",
                  "  - Specify `sample_size` in set_agd_*(), or",
                  "  - Specify covariate values for relative effects using the `newdata` argument",
                  sep = "\n"))

        if (has_agd_arm(object$network)) {
          if (inherits(object$network, "mlnmr_data")) {
            dat_agd_arm <- .unnest_integration(object$network$agd_arm) %>%
              dplyr::mutate(.sample_size = .data$.sample_size / object$network$n_int)
          } else {
            dat_agd_arm <- object$network$agd_arm
          }

          # Only take necessary columns
          dat_agd_arm <- get_model_data_columns(dat_agd_arm, regression = object$regression, label = "AgD (arm-based)")
        } else {
          dat_agd_arm <- tibble::tibble()
        }

        if (has_ipd(object$network)) {
          dat_ipd <- object$network$ipd

          # Only take necessary columns
          dat_ipd <- get_model_data_columns(dat_ipd, regression = object$regression, label = "IPD")

          dat_ipd$.sample_size <- 1
        } else {
          dat_ipd <- tibble::tibble()
        }

        preddat <- dplyr::bind_rows(dat_ipd, dat_agd_arm)
      }

      preddat <- preddat %>%
        dplyr::rename(.trt_old = .data$.trt) %>%
        dplyr::left_join(tidyr::expand(., .study = .data$.study,
                                          .trt = .data$.trt_old),
                         by = ".study")

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
                                      classes = !is.null(object$network$classes))
      X_all <- X_list$X_ipd
      rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

      offset_all <- X_list$offset_ipd

      # Get posterior samples
      post <- as.array(object, pars = c("mu", "d", "beta"))

    # With baseline and newdata specified
    } else {

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

      preddat$.sample_size <- 1

      # Make design matrix of all studies and all treatments
      if (rlang::has_name(preddat, ".trt")) preddat <- dplyr::select(preddat, -.data$.trt)
      preddat <- dplyr::left_join(preddat,
                                  tidyr::expand(preddat,
                                                .study = .data$.study,
                                                .trt = object$network$treatments),
                                  by = ".study")

      # Add in .trtclass if defined in network
      if (!is.null(object$network$classes)) {
        preddat$.trtclass <- object$network$classes[as.numeric(preddat$.trt)]
      }

      # Design matrix, just treating all data as IPD
      X_list <- make_nma_model_matrix(nma_formula,
                                      dat_ipd = preddat,
                                      xbar = object$xbar,
                                      consistency = object$consistency,
                                      classes = !is.null(object$network$classes))
      X_all <- X_list$X_ipd
      rownames(X_all) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

      offset_all <- X_list$offset_ipd

      # Get posterior samples
      post_temp <- as.array(object, pars = c("d", "beta"))

      # Generate baseline samples
      dim_post_temp <- dim(post_temp)
      dim_mu <- c(dim_post_temp[1:2], dplyr::n_distinct(preddat$.study))
      u <- runif(prod(dim_mu))
      mu <- array(rlang::eval_tidy(rlang::call2(baseline$qfun, p = u, !!! baseline$args)),
                  dim = dim_mu)

      # Convert to samples on network ref trt if trt_ref given
      if (!is.null(trt_ref) && trt_ref != nrt) {
        mu <- sweep(mu, 1:2, post_temp[ , , paste0("d[", trt_ref, "]"), drop = FALSE], FUN = "-")
      }

      # Combine mu, d, and beta
      dim_post <- c(dim_post_temp[1:2], dim_mu[3] + dim_post_temp[3])
      post <- array(NA_real_, dim = dim_post)
      post[ , , 1:dim_mu[3]] <- mu
      post[ , , dim_mu[3] + 1:dim_post_temp[3]] <- post_temp

    }

    # Get cutoffs for ordered models
    if (object$likelihood == "ordered") {
      cc <- as.array(object, pars = "cc")
      n_cc <- dim(cc)[3]
      l_cc <- stringr::str_replace(dimnames(cc)[[3]], "^cc\\[(.+)\\]$", "\\1")
    }

    if (level == "individual") {

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

          X_weighted_mean <- matrix(0, ncol = dim(s_pred_array)[3], nrow = n_trt * n_cc)

          X_weighted_mean[cbind(dplyr::group_indices(s_preddat),
                                1:dim(s_pred_array)[3])] <- s_preddat$.weights
        } else {
          s_preddat <- preddat[ss, c(".study", ".trt", ".sample_size")] %>%
            dplyr::mutate(.study = forcats::fct_inorder(forcats::fct_drop(.data$.study)),
                          .trt = forcats::fct_inorder(.data$.trt)) %>%
            dplyr::group_by(.data$.study, .data$.trt) %>%
            dplyr::mutate(.weights = .data$.sample_size / sum(.data$.sample_size))

          X_weighted_mean <- matrix(0, ncol = dim(s_pred_array)[3], nrow = n_trt)

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
                                           .before = 1)
      } else {
        pred_summary <- tibble::add_column(pred_summary, .study = preddat$.study, .before = 1)
      }
      out <- list(summary = pred_summary, sims = pred_array)
    } else {
      out <- list(sims = pred_array)
    }

  }

  if (summary) {
    if (object$likelihood == "ordered")
      class(out) <- c("ordered_nma_summary", "nma_summary")
    else
      class(out) <- "nma_summary"
    attr(out, "xlab") <- "Treatment"
    attr(out, "ylab") <- get_scale_name(likelihood = object$likelihood,
                                        link = object$link,
                                        measure = "absolute",
                                        type = type)
  }
  return(out)
}
