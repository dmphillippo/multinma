#' Relative treatment effects
#'
#' Generate (population-average) relative treatment efects. If a ML-NMR or
#' meta-regression model was fitted, these are specific to each study
#' population.
#'
#' @param x A `stan_nma` object created by [nma()]
#' @param newdata Only used if a regression model is fitted. A data frame of
#'   study details, one row per study, giving the covariate values at which to
#'   produce relative effects. Column names must match variables in the
#'   regression model. If `NULL`, relative effects are produced for all studies
#'   in the network.
#' @param study Column of `newdata` which specifies study names, otherwise
#'   studies will be labelled by row number.
#' @param all_contrasts Logical, generate estimates for all contrasts (`TRUE`),
#'   or just the "basic" contrasts against the network reference treatment
#'   (`FALSE`)? Default `FALSE`.
#'
#' @return A [nma_summary] object.
#' @export
#'
#' @examples
relative_effects <- function(x, newdata = NULL, study = NULL, all_contrasts = FALSE) {

  # Checks
  if (!inherits(x, "stan_nma")) abort("Expecting a `stan_nma` object, as returned by nma().")

  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) abort("`newdata` is not a data frame.")

    if (rlang::quo_is_null(study)) {
      newdata$.study <- 1:nrow(newdata)
    } else {
      newdata <- dplyr::mutate(newdata, .study = {{ study }})

      if (anyDuplicated(newdata$.study))
        abort("Duplicate values in `study` column. Expecting one row per study.")
    }
  }

  if (!is.logical(all_contrasts) || length(all_contrasts) > 1)
    abort("`all_contrasts` should be TRUE or FALSE.")

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce relative effects under inconsistency '{x$consistency}' model."))

  # Produce relative effects
  if (is.null(x$regression)) {
    # If no regression model, relative effects are just the d's

  } else {
    # If regression model, relative effects are study-specific

    if (is.null(newdata)) {
      # Produce relative effects for all studies in network

    } else {
      # Produce relative effects for all studies in newdata

    }
  }
}
