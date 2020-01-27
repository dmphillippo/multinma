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

#' Make all treatment contrasts
#'
#' @param d A 3D MCMC array of basic treatment effects
#' @param trt_ref String naming the reference treatment
#'
#' @return A 3D MCMC array of all contrasts
#' @noRd
make_all_contrasts <- function(d, trt_ref) {
  if (!is.array(d) || length(dim(d)) != 3) abort("Not a 3D MCMC array [Iterations, Chains, Treatments]")
  if (!rlang::is_string(trt_ref)) abort("`trt_ref` must be a single string")

  trts <- c(trt_ref, stringr::str_extract(dimnames(d)[[3]], "(?<=\\[)(.+)(?=\\]$)"))
  ntrt <- length(trts)

  d_ab <- utils::combn(ntrt, 2)

  dim_out <- dim(d)
  dim_out[3] <- ncol(d_ab)

  contrs <- array(NA_real_, dim = dim_out)

  contrs[ , , 1:(ntrt - 1)] <- d
  for (i in ntrt:ncol(d_ab)) {
    contrs[ , , i] <- d[ , , d_ab[2, i] - 1] - d[ , , d_ab[1, i] - 1]
  }

  new_dimnames <- dimnames(d)
  new_dimnames[[3]] <- paste0("d[", trts[d_ab[2, ]], " vs. ", trts[d_ab[1, ]], "]")
  dimnames(contrs) <- new_dimnames

  return(contrs)
}
