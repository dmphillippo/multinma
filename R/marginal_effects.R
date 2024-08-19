#' Marginal treatment effects
#'
#' Generate population-average marginal treatment effects. These are formed from
#' population-average absolute predictions, so this function is a wrapper around
#' [predict.stan_nma()].
#'
#' @param object A `stan_nma` object created by [nma()].
#' @param ... Arguments passed to [predict.stan_nma()], for example to specify
#'   the covariate distribution and baseline risk for a target population, e.g.
#'   `newdata`, `baseline`, and related arguments. For survival outcomes, `type`
#'   can also be specified to determine the quantity from which to form a
#'   marginal effect. For example, `type = "hazard"` with `mtype = "ratio"`
#'   produces marginal hazard ratios, `type = "median"` with `mtype =
#'   "difference"` produces marginal median survival time differences, and so
#'   on.
#' @param mtype The type of marginal effect to construct from the average
#'   absolute effects, either `"difference"` (the default) for a difference of
#'   absolute effects such as a risk difference, `"ratio"` for a ratio of
#'   absolute effects such as a risk ratio, or `"link"` for a difference on the
#'   scale of the link function used in fitting the model such as a marginal log
#'   odds ratio.
#' @param all_contrasts Logical, generate estimates for all contrasts (`TRUE`),
#'   or just the "basic" contrasts against the network reference treatment
#'   (`FALSE`)? Default `FALSE`.
#' @param trt_ref Reference treatment to construct relative effects against, if
#'   `all_contrasts = FALSE`. By default, relative effects will be against the
#'   network reference treatment. Coerced to character string.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param predictive_distribution Logical, when a random effects model has been
#'   fitted, should the predictive distribution for marginal effects in a new
#'   study be returned? Default `FALSE`.
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @examples ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Marginal risk difference in each study population in the network
#' marginal_effects(smk_fit_RE, mtype = "difference")
#'
#' # Since there are no covariates in the model, the marginal and conditional
#' # (log) odds ratios here coincide
#' marginal_effects(smk_fit_RE, mtype = "link")
#' relative_effects(smk_fit_RE)
#'
#' # Marginal risk differences in a population with 67 observed events out of
#' # 566 individuals on No Intervention, corresponding to a Beta(67, 566 - 67)
#' # distribution on the baseline probability of response
#' (smk_rd_RE <- marginal_effects(smk_fit_RE,
#'                                baseline = distr(qbeta, 67, 566 - 67),
#'                                baseline_type = "response",
#'                                mtype = "difference"))
#' plot(smk_rd_RE)
#' }
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_mlnmr_example
#' @examples \donttest{
#' # Population-average marginal probit differences in each study in the network
#' (pso_marg <- marginal_effects(pso_fit, mtype = "link"))
#' plot(pso_marg, ref_line = c(0, 1))
#'
#' # Population-average marginal probit differences in a new target population,
#' # with means and SDs or proportions given by
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
#' # Population-average marginal probit differences of achieving PASI 75 in this
#' # target population, given a Normal(-1.75, 0.08^2) distribution on the
#' # baseline probit-probability of response on Placebo (at the reference levels
#' # of the covariates), are given by
#' (pso_marg_new <- marginal_effects(pso_fit,
#'                                   mtype = "link",
#'                                   newdata = new_agd_int,
#'                                   baseline = distr(qnorm, -1.75, 0.08)))
#' plot(pso_marg_new)
#' }
#'
#' ## Progression free survival with newly-diagnosed multiple myeloma
#' @template ex_ndmm_example
#' @examples \donttest{
#' # We can produce a range of marginal effects from models with survival
#' # outcomes, specified with the mtype and type arguments. For example:
#'
#' # Marginal survival probability difference at 5 years, all contrasts
#' marginal_effects(ndmm_fit, type = "survival", mtype = "difference",
#'                  times = 5, all_contrasts = TRUE)
#'
#' # Marginal difference in RMST up to 5 years
#' marginal_effects(ndmm_fit, type = "rmst", mtype = "difference", times = 5)
#'
#' # Marginal median survival time ratios
#' marginal_effects(ndmm_fit, type = "median", mtype = "ratio")
#'
#' # Marginal log hazard ratios
#' # With no covariates in the model, these are constant over time and study
#' # populations, and are equal to the log hazard ratios from relative_effects()
#' plot(marginal_effects(ndmm_fit, type = "hazard", mtype = "link"),
#'      # The hazard is infinite at t=0 in some studies, giving undefined logHRs at t=0
#'      na.rm = TRUE)
#'
#' # The NDMM vignette demonstrates the production of time-varying marginal
#' # hazard ratios from a ML-NMR model that includes covariates, see
#' # `vignette("example_ndmm")`
#'
#' # Marginal survival difference over time
#' plot(marginal_effects(ndmm_fit, type = "survival", mtype = "difference"))
#' }
#'
marginal_effects <- function(object,
                             ...,
                             mtype = c("difference", "ratio", "link"),
                             all_contrasts = FALSE, trt_ref = NULL,
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                             predictive_distribution = FALSE,
                             summary = TRUE) {

  # Checks
  if (!inherits(object, "stan_nma")) abort("Expecting a `stan_nma` object, as returned by nma().")

  mtype <- rlang::arg_match(mtype)

  if (!rlang::is_bool(all_contrasts))
    abort("`all_contrasts` should be TRUE or FALSE.")

  if (!is.null(trt_ref)) {
    if (all_contrasts) {
      warn("Ignoring `trt_ref` when all_contrasts = TRUE.")
      trt_ref <- levels(object$network$treatments)[1]
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
    trt_ref <- levels(object$network$treatments)[1]
  }

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  check_probs(probs)

  if(!rlang::is_bool(predictive_distribution))
    abort("`predictive_distribution` should be TRUE or FALSE")
  if (predictive_distribution && object$trt_effects != "random") predictive_distribution <- FALSE

  # Cannot produce marginal effects for inconsistency models
  if (object$consistency != "consistency")
    abort(glue::glue("Cannot produce marginal effects under inconsistency '{object$consistency}' model."))

  # Call predict to get absolute predictions
  if ("level" %in% ...names()) warn('Ignoring `level` argument, this is always "aggregate" for marginal effects.')

  if (object$likelihood %in% valid_lhood$survival) {
    pred <- rlang::eval_tidy(rlang::call2(
              stats::predict,
              !!! rlang::dots_list(object,
                                   level = "aggregate",
                                   summary = FALSE,
                                   predictive_distribution = predictive_distribution,
                                   !!! rlang::enquos(...),
                                   .homonyms = "first")))
  } else {
    if ("type" %in% ...names()) warn('Ignoring `type` argument, this is always "response" for marginal effects (except for survival models).')
    pred <- rlang::eval_tidy(rlang::call2(
              stats::predict,
              !!! rlang::dots_list(object,
                                   type = "response", level = "aggregate",
                                   summary = FALSE,
                                   predictive_distribution = predictive_distribution,
                                   !!! rlang::enquos(...),
                                   .homonyms = "first")))
  }

  # Set up output metadata
  pred_meta <- pred$summary
  vars <- intersect(c(".study", ".time", ".quantile", ".category"), colnames(pred_meta))

  if (".time" %in% vars) {
    pred_meta <- { if (".study" %in% vars) dplyr::group_by(pred_meta, .data$.study, .data$.trt) else dplyr::group_by(pred_meta, .data$.trt) } %>%
      dplyr::mutate(.time_id = 1:dplyr::n()) %>%
      dplyr::ungroup()

    vars <- c(vars, ".time_id")
  }

  out_meta <- pred_meta

  if (!all_contrasts) {
    out_meta <- dplyr::filter(out_meta, .data$.trt != trt_ref) %>%
      dplyr::mutate(.trtb = .data$.trt, .trta = factor(trt_ref, levels = levels(object$network$treatments))) %>%
      dplyr::select(dplyr::any_of(".study"), ".trtb", ".trta", dplyr::everything(), -".trt")
  } else {
    contrs <- utils::combn(object$network$treatments, 2)
    trtb <- contrs[2, ]
    trta <- contrs[1, ]
    out_meta <- dplyr::distinct(out_meta, dplyr::pick(dplyr::all_of(vars))) %>%
      dplyr::mutate(.trtb = list(trtb), .trta = list(trta))  %>%
      tidyr::unnest(cols = c(".trtb", ".trta")) %>%
      dplyr::select(dplyr::any_of(".study"), ".trtb", ".trta", dplyr::everything())
  }

  pred_meta$id <- 1:nrow(pred_meta)
  out_meta$id <- 1:nrow(out_meta)

  # Output parameter names
  pnames <- paste0("marg[",
                   if (rlang::has_name(out_meta, ".study")) paste0(out_meta$.study, ": ") else character(),
                   out_meta$.trtb,
                   if (all_contrasts) paste0(" vs. ", out_meta$.trta) else character(),
                   if (".time" %in% vars) paste0(", ", out_meta$.time_id)
                   else if (".quantile" %in% vars) paste0(", ", out_meta$.quantile)
                   else if (".category" %in% vars) paste0(", ", out_meta$.category)
                   else character()
                   , "]")

  # Set up output array
  d_out <- dim(pred$sims)
  d_out[[3]] <- nrow(out_meta)
  dn_out <- list(iterations = NULL, chains = NULL, parameters = pnames)
  out_array <- array(NA_real_, dim = d_out, dimnames = dn_out)

  # Set contrast function
  if (mtype == "difference") {
    cf <- function(b, a) b - a
  } else if (mtype == "ratio") {
    cf <- function(b, a) b / a
  } else if (mtype == "link") {
    cf <- function(b, a) link_fun(b, link = object$link) - link_fun(a, link = object$link)
  }

  cfv <- function(x) cf(x[2], x[1])
  mk_contr <- function(x) utils::combn(x, 2, FUN = cfv)

  # Calculate contrasts, looping over groups
  grps <- dplyr::distinct(out_meta, dplyr::pick(dplyr::all_of(vars)))

  for (i in 1:ifelse(length(vars), nrow(grps), 1)) {
    if (length(vars)) {
      grpi <- grps[i, ]
      predi <- dplyr::semi_join(pred_meta, grpi, by = vars)$id
      outi <- dplyr::semi_join(out_meta, grpi, by = vars)$id
    } else { # No grouping variables
      grpi <- tibble::tibble(.rows = 1)
      predi <- pred_meta$id
      outi <- out_meta$id
    }

    if (all_contrasts && nlevels(object$network$treatments) > 2) {
      out_array[ , , outi] <- aperm(apply(pred$sims[ , , predi], MARGIN = 1:2, FUN = mk_contr), c(2, 3, 1))
    } else {
      predi_ref <- dplyr::semi_join(pred_meta, dplyr::mutate(grpi, .trt = trt_ref), by = c(vars, ".trt"))$id
      predi_nonref <- setdiff(predi, predi_ref)
      out_array[ , , outi] <- sweep(pred$sims[ , , predi_nonref, drop = FALSE], MARGIN = 1:2, STATS = pred$sims[ , , predi_ref], FUN = cf)
    }
  }

  out_meta <- dplyr::select(out_meta, -"id")
  if (".time_id" %in% vars) out_meta <- dplyr::select(out_meta, -".time_id")

  # Summarise
  if (summary) {
    out_summary <- summary_mcmc_array(out_array, probs)
    out_summary <- dplyr::bind_cols(out_meta, out_summary)
  } else {
    out_summary <- out_meta
  }

  out <- list(summary = out_summary, sims = out_array)

  # Set attributes
  if (summary) {
    attr(out, "xlab") <- "Treatment"

    if (object$likelihood %in% valid_lhood$survival) {
      dots <- rlang::enquos(...)
      if (rlang::has_name(dots, "type")) type <- rlang::eval_tidy(dots$type)
      else type <- "survival"
    } else if (mtype == "link") {
      type <- "link"
    } else {
      type <- "response"
    }

    if (mtype == "difference") {
      attr(out, "ylab") <- paste0("Marginal ",
                                  get_scale_name(likelihood = object$likelihood,
                                                 link = object$link,
                                                 measure = "absolute",
                                                 type = type),
                                  " Difference")
    } else if (mtype == "ratio") {
      attr(out, "ylab") <- paste0("Marginal ",
                                  get_scale_name(likelihood = object$likelihood,
                                                 link = object$link,
                                                 measure = "absolute",
                                                 type = type),
                                  " Ratio")
    } else if (mtype == "link") {
      if (object$likelihood %in% valid_lhood$survival) {
        if (object$link == "log") {
          attr(out, "ylab") <- paste0("Marginal log ",
                                      get_scale_name(likelihood = object$likelihood,
                                                     link = object$link,
                                                     measure = "absolute",
                                                     type = type),
                                      " Ratio")
        } else {
          attr(out, "ylab") <- paste0("Marginal ", object$link, " ",
                                      get_scale_name(likelihood = object$likelihood,
                                                     link = object$link,
                                                     measure = "absolute",
                                                     type = type),
                                      " Difference ")
        }
      } else {
        attr(out, "ylab") <- paste0("Marginal ",
                                    get_scale_name(likelihood = object$likelihood,
                                                   link = object$link,
                                                   measure = "relative",
                                                   type = type))
      }
    }

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
