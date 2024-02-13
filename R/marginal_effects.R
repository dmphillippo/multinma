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
  }

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  check_probs(probs)

  if(!rlang::is_bool(predictive_distribution))
    abort("`predictive_distribution` should be TRUE or FALSE")
  if (predictive_distribution && x$trt_effects != "random") predictive_distribution <- FALSE

  # Cannot produce marginal effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce marginal effects under inconsistency '{x$consistency}' model."))

  # Get network reference treatment
  nrt <- levels(x$network$treatments)[1]

  # Call predict to get absolute predictions
  pred <- predict(object,
                  type = "response", level = "aggregate",
                  summary = FALSE,
                  predictive_distribution = predictive_distribution,
                  ...)
}
