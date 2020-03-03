#' Predictions of absolute effects from NMA models
#'
#' Obtain predictions of absolute effects from NMA models fitted with [nma()].
#' For example, if a model is fitted to binary data with a logit link, predicted
#' outcome probabilities or log odds can be produced.
#'
#' @param object A `stan_nma` object created by [nma()]
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
#'   If `NULL`, prections are produced for all studies with IPD and/or
#'   arm-based AgD in the network, depending on the value of `type`.
#' @param study Column of `newdata` which specifies study names or IDs. When not
#'   specified: if `newdata` contains integration points produced by
#'   [add_integration()], studies will be labelled sequentially by row;
#'   otherwise data will be assumed to come from a single study.
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
#' @examples
predict.stan_nma <- function(object,
                             baseline = NULL, newdata = NULL, study = NULL,
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
      newdata <- dplyr::mutate(newdata, .study = nfactor(.study))
    }
  }

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  # Cannot produce predictions for inconsistency models
  if (object$consistency != "consistency")
    abort(glue::glue("Cannot produce predictions under inconsistency '{x$consistency}' model."))

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
        studies <- forcats::fct_unique(dplyr::bind_rows(object$network$ipd, object$network$agd_arm)$.study)
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

      # Combine mu and d
      dim_post <- c(dim_d[1:2], dim_d[3] + 1)
      post <- array(NA_real_, dim = dim_post)
      post[ , , 1] <- mu
      post[ , , 2:dim_post[3]] <- d

      # Get prediction array
      pred_array <- tcrossprod_mcmc_array(post, X_all)

    }

    # Transform predictions if type = "response"
    if (type == "response") {
      pred_array <- inverse_link(pred_array, link = object$link)
    }

    # Produce nma_summary
    if (summary) {
      pred_summary <- summary_mcmc_array(pred_array, probs)
      if (is.null(baseline))
        pred_summary <- tibble::add_column(pred_summary, .study = preddat$.study, .before = 1)
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
        } else {
          dat_agd_arm <- tibble::tibble()
        }

        if (has_ipd(object$network)) {
          dat_ipd <- object$network$ipd
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
      post <- as.array(object, pars = c("mu", "d", "beta"))

    # With baseline and newdata specified
    } else {

      if (level == "individual") {
        if (!has_ipd(object$network))
          warn("Producing individual predictions from an aggregate-level regression. Interpret with great caution!")

        preddat <- newdata
      } else {
        if (inherits(object, "stan_mlnmr")) {
          if (!inherits(newdata, "integration_df")) {
            abort("No integration points found in `newdata`. Specify integration points using add_integration().")
          } else {
            preddat <- .unnest_integration(newdata)
          }
        } else {
          if (has_ipd(object$network) && inherits(newdata, "integration_df")) {
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

      # Combine mu, d, and beta
      dim_post <- c(dim_post_temp[1:2], dim_mu[3] + dim_post_temp[3])
      post <- array(NA_real_, dim = dim_post)
      post[ , , 1:dim_mu[3]] <- mu
      post[ , , dim_mu[3] + 1:dim_post_temp[3]] <- post_temp

    }

    # Get prediction array
    pred_array <- tcrossprod_mcmc_array(post, X_all)

    if (!is.null(offset_all))
      pred_array <- sweep(pred_array, 3, offset_all, FUN = "+")

    # Transform predictions if type = "response"
    if (type == "response") {
      pred_array <- inverse_link(pred_array, link = object$link)
    }

    # Aggregate predictions if level = "aggregate"
    if (level == "aggregate") {
      studies <- preddat$.study
      n_studies <- nlevels(studies)
      n_trt <- nlevels(object$network$treatments)

      preddat <- preddat %>%
        dplyr::group_by(.data$.study, .data$.trt) %>%
        dplyr::mutate(.weights = .data$.sample_size / sum(.data$.sample_size))

      X_weighted_mean <- matrix(0, ncol = dim(pred_array)[3], nrow = n_studies * n_trt)

      X_weighted_mean[cbind(dplyr::group_indices(preddat),
                            1:dim(pred_array)[3])] <- preddat$.weights

      preddat <- dplyr::summarise(preddat, .weights = list(.data$.weights))

      rownames(X_weighted_mean) <- paste0("pred[", preddat$.study, ": ", preddat$.trt, "]")

      pred_array <- tcrossprod_mcmc_array(pred_array, X_weighted_mean)
    }

    # Produce nma_summary
    if (summary) {
      pred_summary <- summary_mcmc_array(pred_array, probs)
      pred_summary <- tibble::add_column(pred_summary, .study = preddat$.study, .before = 1)
      out <- list(summary = pred_summary, sims = pred_array)
    } else {
      out <- list(sims = pred_array)
    }

  }

  if (summary) {
    class(out) <- "nma_summary"
    attr(out, "xlab") <- "Treatment"
    attr(out, "ylab") <- get_scale_name(likelihood = object$likelihood,
                                        link = object$link,
                                        measure = "absolute",
                                        type = type)
  }
  return(out)
}
