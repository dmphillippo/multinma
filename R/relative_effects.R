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

    .study <- pull_non_null(newdata, enquo(study))
    if (is.null(.study)) {
      newdata$.study <- nfactor(seq_len(nrow(newdata)))
    } else {
      newdata <- dplyr::mutate(newdata, .study = nfactor(.study))

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
    out <- list(summary = re_summary, sims = re_array)

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

      dat_all <- dplyr::bind_rows(dat_agd_arm, dat_agd_contrast, dat_ipd)

      # Take the first row for each study. We will correct the design matrix below
      dat_studies <- dat_all %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::slice(1)

    } else {
      # Produce relative effects for all studies in newdata

      dat_studies <- newdata
    }

    # Get number of treatments
    ntrt <- nlevels(x$network$treatments)

    # Apply centering if used
    if (!is.null(x$xbar)) {
      cen_vars <- intersect(names(dat_studies), names(x$xbar))
      dat_studies[, cen_vars] <- sweep(dat_studies[, cen_vars, drop = FALSE], 2, x$xbar[cen_vars])
    }

    # Expand rows for every treatment
    all_trts <- tidyr::expand_grid(.study = dat_studies$.study, .trt = x$network$treatments[-1])
    if (rlang::has_name(dat_studies, ".trt")) dat_studies <- dplyr::select(dat_studies, -.data$.trt)
    dat_studies <- dplyr::left_join(all_trts, dat_studies, by = ".study")

    # Add in .trtclass if defined in network
    if (!is.null(x$network$classes)) {
      dat_studies$.trtclass <- x$network$classes[as.numeric(dat_studies$.trt)]
    }

    # Sanitise factor levels
    dat_studies <- dplyr::mutate_at(dat_studies,
      .vars = if (!is.null(x$network$classes)) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
      .funs = fct_sanitise
    )

    # Get model formula and design matrix
    nma_formula <- make_nma_formula(x$regression,
                                    consistency = x$consistency,
                                    classes = !is.null(x$network$classes),
                                    class_interactions = x$class_interactions)


    # Drop study to factor to 1L if only one study (avoid contrasts need 2 or
    # more levels error)
    if (dplyr::n_distinct(dat_studies$.study) == 1) {

      # Save study label to restore
      single_study_label <- unique(dat_studies$.study)
      dat_studies$.study_temp <- dat_studies$.study
      dat_studies$.study <- 1L

      # Fix up model formula with an intercept
      nma_formula <- update.formula(nma_formula, ~. + 1)
    } else {
      single_study_label <- NULL
    }

    X_all <- model.matrix(nma_formula, data = dat_studies)

    if (!is.null(single_study_label)) {
      # Restore single study label and .study column
      colnames(X_all) <- stringr::str_replace(colnames(X_all),
                                              "^\\.study$",
                                              paste0(".study", single_study_label))
      dat_studies <- dat_studies %>%
        dplyr::mutate(.study = .data$.study_temp) %>%
        dplyr::select(-.data$.study_temp)

      # Drop intercept column from design matrix
      X_all <- X_all[, -1, drop = FALSE]
    }

    # Remove columns for reference level of .trtclass
    if (!is.null(x$network$classes)) {
      col_trtclass_ref <- grepl(paste0(".trtclass", levels(x$network$classes)[1]),
                                colnames(X_all), fixed = TRUE)
      X_all <- X_all[, !col_trtclass_ref, drop = FALSE]
    }

    # Subset design matrix into EM columns and trt columns
    X_d <- X_all[, grepl("^(\\.trt|\\.contr)[^:]+$", colnames(X_all))]
    EM_regex <- "(^\\.trt(class)?.+\\:)|(\\:\\.trt(class)?.+$)"
    X_EM <- X_all[, grepl(EM_regex, colnames(X_all))]

    # If there are no EMs (regression model had no interactions with .trt) then
    # just return the treatment effects (not study specific)
    if (ncol(X_EM) == 0) {
      re_array <- as.array(as.stanfit(x), pars = "d")
      if (all_contrasts) {
        re_array <- make_all_contrasts(re_array, trt_ref = trt_ref)
      }
      re_summary <- summary_mcmc_array(re_array, probs = probs)
      out <- list(summary = re_summary, sims = re_array)
    } else {

      # Figure out which covariates are EMs from model matrix
      EM_col_names <- stringr::str_remove(colnames(X_EM), EM_regex)
      EM_vars <- unique(EM_col_names)

      # Replace EM design matrix with study means if newdata is NULL
      if (is.null(newdata)) {

        # Apply centering if used
        if (!is.null(x$xbar)) {
          cen_vars <- intersect(names(dat_all), names(x$xbar))
          dat_all_cen <- dat_all
          dat_all_cen[, cen_vars] <- sweep(dat_all[, cen_vars, drop = FALSE], 2, x$xbar[cen_vars])
        } else {
          dat_all_cen <- dat_all
        }

        # Get model matrix of EM "main effects" - notably this expands out factors
        # into dummy variables so we can average those too
        EM_formula <- as.formula(paste0("~", paste(EM_vars, collapse = " + ")))


        # Calculate mean covariate values by study in the network
        X_study_means <- model.matrix(EM_formula, data = dat_all_cen) %>%
          tibble::as_tibble() %>%
          tibble::add_column(.study = dat_all$.study, .sample_size = dat_all$.sample_size, .before = 1) %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::summarise_at(setdiff(colnames(.), c(".sample_size", ".study")),
                              ~weighted.mean(., w = .data$.sample_size)) %>%
          dplyr::select(!!! unique(EM_col_names)) %>%
          as.matrix()

        # Repeat columns across interaction columns in X_EM
        X_study_means_rep <- X_study_means[, EM_col_names, drop = FALSE]

        # Repeat study rows for the number of treatment parameters
        X_study_means_rep <- apply(X_study_means_rep, 2, rep, each = ntrt - 1)

        # Replace non-zero entries of design matrix X_EM with corresponding mean values
        # This works only because trt columns are 0/1, so interactions are just the covariate values
        nonzero <- X_EM != 0
        X_EM[nonzero] <- X_study_means_rep[nonzero]
      }

      # Name columns to match Stan parameters
      colnames(X_EM) <- paste0("beta[", colnames(X_EM), "]")
      colnames(X_d) <- paste0("d[", x$network$treatments[-1], "]")

      X_EM_d <- cbind(X_EM, X_d)

      # Name rows by treatment for now (required for make_all_contrasts)
      rownames(X_EM_d) <- paste0("d[", rep(x$network$treatments[-1],
                                           times = dplyr::n_distinct(dat_studies$.study)) , "]")

      # Linear combination with posterior MCMC array
      d_array <- as.array(x, pars = colnames(X_EM_d))
      re_array <- tcrossprod_mcmc_array(d_array, X_EM_d)

      # Produce all contrasts, if required
      if (all_contrasts)
        re_array <- make_all_contrasts(re_array, trt_ref = levels(x$network$treatments)[1])

      # Add in study names to parameters
      parnames <- stringr::str_extract(dimnames(re_array)[[3]], "(?<=^d\\[)(.+)(?=\\]$)")
      study_parnames <- paste0("d[", dat_studies$.study, ": ", parnames, "]")
      dimnames(re_array)[[3]] <- study_parnames

      # Create summary stats
      re_summary <- summary_mcmc_array(re_array, probs = probs) %>%
        tibble::add_column(.study =
                             if (all_contrasts) {
                               rep(unique(dat_studies$.study), each = ntrt * (ntrt - 1) / 2)
                             } else {
                               dat_studies$.study
                             },
                           .before = 1)

      # Prepare study covariate info
      if (is.null(newdata)) {
        study_EMs <- X_study_means

        # Uncenter if necessary
        if (!is.null(x$xbar)) {
          cen_vars <- intersect(colnames(study_EMs), names(x$xbar))
          study_EMs[, cen_vars] <- sweep(study_EMs[, cen_vars, drop = FALSE], 2, x$xbar[cen_vars], FUN = "+")
        }
      } else {
        study_EMs <- newdata[EM_vars]
      }

      study_EMs <- tibble::as_tibble(study_EMs) %>%
        tibble::add_column(.study = unique(dat_studies$.study), .before = 1)

      out <- list(summary = re_summary, sims = re_array, studies = study_EMs)
    }
  }

  # Return nma_summary object
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
#' @noRd
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
