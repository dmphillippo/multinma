#' Marginal treatment effects
#'
#' Generate population-average marginal treatment effects. These are formed from
#' population-average absolute predictions, so this function is a wrapper around
#' [predict.stan_nma()].
#'
#' @param object A `stan_nma` object created by [nma()].
#' @param ... Arguments passed to [predict.stan_nma()], for example to specify
#'   the covariate distribution and baseline risk for a target population.
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
#' @return
#' @export
#'
#' @examples
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
  if (object$likelihood %in% valid_lhood$survival) {
    pred <- predict(object,
                    level = "aggregate",
                    summary = FALSE,
                    predictive_distribution = predictive_distribution,
                    ...)
  } else {
    pred <- predict(object,
                    type = "response", level = "aggregate",
                    summary = FALSE,
                    predictive_distribution = predictive_distribution,
                    ...)
  }

  # Set up output metadata
  pred_meta <- pred$summary
  out_meta <- pred_meta

  vars <- intersect(c(".study", ".time", ".quantile", ".category"), colnames(pred_meta))

  if (!all_contrasts) {
    out_meta <- dplyr::filter(out_meta, .data$.trt != trt_ref) %>%
      dplyr::mutate(.trtb = .data$.trt, .trta = trt_ref) %>%
      dplyr::select(-".trt", ".study", ".trtb", ".trta", dplyr::everything())
  } else {
    contrs <- utils::combn(object$network$treatments, 2)
    trtb <- contrs[2, ]
    trta <- contrs[1, ]
    out_meta <- dplyr::distinct(out_meta, dplyr::pick(dplyr::all_of(vars))) %>%
      dplyr::mutate(.trtb = list(trtb), .trta = list(trta))  %>%
      tidyr::unnest(cols = c(.data$.trtb, .data$.trta)) %>%
      dplyr::select(".study", ".trtb", ".trta", dplyr::everything()) %>%
      dplyr::arrange(".study", ".trta", ".trtb", dplyr::pick(dplyr::everything()))
  }

  pred_meta$id <- 1:nrow(pred_meta)
  out_meta$id <- 1:nrow(out_meta)

  # Output parameter names
  pnames <- paste0("D[",
                   if (dplyr::n_distinct(out_meta$.study) > 1) paste0(out_meta$.study, ": ") else character(),
                   out_meta$.trtb,
                   if (all_contrasts) paste0(" vs. ", out_meta$.trta) else character(),
                   if (".time" %in% vars) paste0(", ", rep(1:dplyr::n_distinct(out_meta$.time), times = dplyr::n_distinct(out_meta$.study, out_meta$.trtb, out_meta$.trta)))
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
  mk_contr <- function(x) combn(x, 2, FUN = cfv)

  # Calculate contrasts, looping over groups
  grps <- dplyr::distinct(out_meta, dplyr::pick(dplyr::all_of(vars)))

  for (i in 1:nrow(grps)) {
    grpi <- grps[i, ]
    predi <- dplyr::semi_join(pred_meta, grpi, by = vars)$id
    outi <- dplyr::semi_join(out_meta, grpi, by = vars)$id

    if (all_contrasts) {
      out_array[ , , outi] <- aperm(apply(pred$sims[ , , predi], MARGIN = 1:2, FUN = mk_contr), c(2, 3, 1))
    } else {
      predi_ref <- dplyr::semi_join(pred_meta, dplyr::mutate(grpi, .trt = trt_ref), by = c(vars, ".trt"))$id
      predi_nonref <- setdiff(predi, predi_ref)
      out_array[ , , outi] <- sweep(pred$sims[ , , predi_nonref, drop = FALSE], MARGIN = 1:2, STATS = pred$sims[ , , predi_ref], FUN = cf)
    }
  }

  out_meta <- dplyr::select(out_meta, -"id")

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
      attr(out, "ylab") <- paste0("Marginal ",
                                  get_scale_name(likelihood = object$likelihood,
                                                 link = object$link,
                                                 measure = "relative",
                                                 type = type))
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
