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
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#'
#' @return A [nma_summary] object.
#' @export
#'
#' @examples
relative_effects <- function(x, newdata = NULL, study = NULL, all_contrasts = FALSE,
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {

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

  # Get reference treatment
  trt_ref <- levels(x$network$treatments)[1]

  # Produce relative effects
  if (is.null(x$regression)) {
    # If no regression model, relative effects are just the d's

    re_array <- as.array(as.stanfit(x), pars = "d")

    if (all_contrasts) {
      re_array <- make_all_contrasts(re_array, trt_ref = trt_ref)
    }

    re_summary <- summary_mcmc_array(re_array, probs = probs)

  } else {
    # If regression model, relative effects are study-specific

    if (is.null(newdata)) {
      # Produce relative effects for all studies in network

      # Make data frame of study covariate means
      if ((has_agd_arm(x$network) || has_agd_contrast(x$network)) && !has_agd_sample_size(x$network))
        abort(
          paste("AgD study sample sizes not specified in network, cannot calculate mean covariate values.",
                "  - Specify `sample_size` in set_agd_*(), or",
                "  - Specify covariate values for relative effects using the `newdata` argument",
                sep = "\n"))


      if (has_agd_arm(x$network)) {
        if (inherits(x$network, "mlnmr_data")) {
          dat_agd_arm <- unnest_integration_points(x$network$agd_arm, x$network$int_names) %>%
            dplyr::mutate(.sample_size = .data$.sample_size / x$network$n_int)
        } else {
          dat_agd_arm <- x$network$agd_arm
        }
      } else {
        dat_agd_arm <- tibble::tibble()
      }

      if (has_agd_contrast(x$network)) {
        if (inherits(x$network, "mlnmr_data")) {
          dat_agd_contrast <- unnest_integration_points(x$network$agd_contrast, x$network$int_names) %>%
            dplyr::mutate(.sample_size = .data$.sample_size / x$network$n_int)
        } else {
          dat_agd_contrast <- x$network$agd_contrast
        }
      } else {
        dat_agd_contrast <- tibble::tibble()
      }

      if (has_ipd(x$network)) {
        if (inherits(x$network, "mlnmr_data")) {
          dat_ipd <- unnest_integration_points(x$network$ipd, x$network$int_names)
        } else {
          dat_ipd <- x$network$ipd
        }
        dat_ipd$.sample_size <- 1
      } else {
        dat_ipd <- tibble::tibble()
      }

      dat_all <- dplyr::bind_rows(dat_agd_arm, dat_agd_contrast, dat_ipd) %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::summarise_if(is.numeric, ~weighted.mean(., w = .data$.sample_size))

    } else {
      # Produce relative effects for all studies in newdata

      dat_all <- newdata
    } else {
      # Produce relative effects for all studies in newdata

    }
  }

  out <- list(summary = re_summary, sims = re_array)
  class(out) <- "nma_summary"

  return(out)
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

#' Crossproduct for 3D MCMC arrays
#'
#' @param x A matrix
#' @param a A 3D MCMC array
#'
#' @return A 3D MCMC array
tcrossprod_mcmc_array <- function(a, x) {
  if (!is.array(a) || length(dim(a)) != 3) abort("Not a 3D MCMC array [Iterations, Chains, Parameters]")

  dim_out <- c(dim(a)[1:2], NROW(x))
  dimnames_out <- dimnames(a)
  dimnames_out[[3]] <- if (!is.null(rownames(x))) rownames(x) else paste0("V", 1:NROW(x))
  out <- array(NA_real_, dim = dim_out)

  nchains <- dim(a)[2]

  for (i in 1:nchains) {
    out[ , i, ] <- tcrossprod(a[ , i, ], x)
  }

  dimnames(out) <- dimnames_out

  return(out)
}
