#' Relative treatment effects
#'
#' Generate (population-average) relative treatment effects. If a ML-NMR or
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
#' @param trt_ref Reference treatment to construct relative effects against, if
#'   `all_contrasts = FALSE`. By default, relative effects will be against the
#'   network reference treatment. Coerced to character string.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param predictive_distribution Logical, when a random effects model has been
#'   fitted, should the predictive distribution for relative effects in a new
#'   study be returned? Default `FALSE`.
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @seealso [plot.nma_summary()] for plotting the relative effects.
#'
#' @examples
#' ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Produce relative effects
#' smk_releff_RE <- relative_effects(smk_fit_RE)
#' smk_releff_RE
#' plot(smk_releff_RE, ref_line = 0)
#'
#' # Relative effects for all pairwise comparisons
#' relative_effects(smk_fit_RE, all_contrasts = TRUE)
#'
#' # Relative effects against a different reference treatment
#' relative_effects(smk_fit_RE, trt_ref = "Self-help")
#'
#' # Transforming to odds ratios
#' # We work with the array of relative effects samples
#' LOR_array <- as.array(smk_releff_RE)
#' OR_array <- exp(LOR_array)
#'
#' # mcmc_array objects can be summarised to produce a nma_summary object
#' smk_OR_RE <- summary(OR_array)
#'
#' # This can then be printed or plotted
#' smk_OR_RE
#' plot(smk_OR_RE, ref_line = 1)
#' }
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_mlnmr_example
#' @examples \donttest{
#' # Produce population-adjusted relative effects for all study populations in
#' # the network
#' pso_releff <- relative_effects(pso_fit)
#' pso_releff
#' plot(pso_releff, ref_line = 0)
#'
#' # Produce population-adjusted relative effects for a different target
#' # population
#' new_agd_means <- data.frame(
#'   bsa = 0.6,
#'   prevsys = 0.1,
#'   psa = 0.2,
#'   weight = 10,
#'   durnpso = 3)
#'
#' relative_effects(pso_fit, newdata = new_agd_means)
#' }
relative_effects <- function(x, newdata = NULL, study = NULL,
                             all_contrasts = FALSE, trt_ref = NULL,
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                             predictive_distribution = FALSE,
                             summary = TRUE) {

  # Checks
  if (!inherits(x, "stan_nma")) abort("Expecting a `stan_nma` object, as returned by nma().")

  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) abort("`newdata` is not a data frame.")

    .study <- pull_non_null(newdata, enquo(study))
    if (is.null(.study)) {
      newdata$.study <- nfactor(paste("New", seq_len(nrow(newdata))))
    } else {
      check_study(.study)
      newdata <- dplyr::mutate(newdata, .study = nfactor(.study))

      if (anyDuplicated(newdata$.study))
        abort("Duplicate values in `study` column. Expecting one row per study.")
    }
  }

  if (!rlang::is_bool(all_contrasts))
    abort("`all_contrasts` should be TRUE or FALSE.")

  if (!is.null(trt_ref)) {
    if (all_contrasts) {
      warn("Ignoring `trt_ref` when all_contrasts = TRUE.")
      trt_ref <- NULL
    } else {
      if (length(trt_ref) > 1) abort("`trt_ref` must be length 1.")
      trt_ref <- as.character(trt_ref)
      lvls_trt <- levels(x$network$treatments)
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

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce relative effects under inconsistency '{x$consistency}' model."))

  # Get network reference treatment
  nrt <- levels(x$network$treatments)[1]

  # Produce relative effects
  if (is.null(x$regression) || is_only_offset(x$regression)) {
    # If no regression model, relative effects are just the d's

    if (!predictive_distribution) {
      re_array <- as.array(as.stanfit(x), pars = "d")
    } else {
      # For predictive distribution, use delta_new instead of d
      re_array <- get_delta_new(x)
    }

    # Set treatments vector
    trtb <- x$network$treatments[-1]
    trta <- x$network$treatments[1]

    if (all_contrasts) {
      re_array <- make_all_contrasts(re_array, trt_ref = nrt)

      contrs <- utils::combn(x$network$treatments, 2)
      trtb <- contrs[2, ]
      trta <- contrs[1, ]

    } else if (!is.null(trt_ref) && trt_ref != nrt) {
      d_ref <- re_array[ , , paste0("d[", trt_ref, "]"), drop = FALSE]
      re_array <- sweep(re_array, 1:2, d_ref, FUN = "-")

      # Add in parameter for network ref trt in place of trt_ref
      re_array[ , , paste0("d[", trt_ref, "]")] <- -d_ref
      d_names <- dimnames(re_array)[[3]]
      d_names[d_names == paste0("d[", trt_ref, "]")] <- paste0("d[", nrt, "]")
      dimnames(re_array)[[3]] <- d_names

      # Reorder parameters
      d_names <- c(paste0("d[", nrt, "]"), d_names[d_names != paste0("d[", nrt, "]")])
      re_array <- re_array[ , , d_names, drop = FALSE]

      trt_reorder <- sort(forcats::fct_relevel(x$network$treatments, trt_ref))
      trtb <- trt_reorder[-1]
      trta <- trt_reorder[1]
    }

    # Fix up parameter names for predictive_distribution = TRUE
    if (predictive_distribution) {
      dimnames(re_array)[[3]] <- gsub("^d\\[", "delta_new[", dimnames(re_array)[[3]])
    }

    if (summary) {
      re_summary <- summary_mcmc_array(re_array, probs = probs) %>%
        tibble::add_column(.trtb = trtb, .trta = trta, .before = 1)

      out <- list(summary = re_summary, sims = re_array)
    } else {
      out <- list(sims = re_array)
    }

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
        if (x$network$outcome$agd_arm != "survival") {
          if (inherits(x$network, "mlnmr_data")) {
            dat_agd_arm <- .unnest_integration(x$network$agd_arm) %>%
              dplyr::mutate(.sample_size = .data$.sample_size / x$network$n_int)
          } else {
            dat_agd_arm <- x$network$agd_arm
          }
        } else {
          # For survival outcomes, check if we need to access .data_orig
          if (any(rlang::has_name(x$agd_arm$.data_orig, all.vars(x$regression)))) {
            dat_agd_arm <- x$network$agd_arm %>%
              # Drop duplicated names in outer dataset from .data_orig before unnesting
              dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(x$network$agd_arm)))),
                            # Reset sample size for weighted mean later
                            .sample_size = 1) %>%
              # Unnest .data_orig
              tidyr::unnest(cols = c(".Surv", ".data_orig"))
          } else {
            # Drop nested .Surv column, breaks join with ipd
            dat_agd_arm <- dplyr::select(x$network$agd_arm, -.data$.Surv)
          }

          # Unnest integration points if present
          if (inherits(x$network, "mlnmr_data")) {
            dat_agd_arm <- .unnest_integration(dat_agd_arm) %>%
              dplyr::mutate(.sample_size = .data$.sample_size / x$network$n_int)
          }
        }

        # Only take necessary columns
        dat_agd_arm <- get_model_data_columns(dat_agd_arm, regression = x$regression, label = "AgD (arm-based)")
      } else {
        dat_agd_arm <- tibble::tibble()
      }

      if (has_agd_contrast(x$network)) {
        if (inherits(x$network, "mlnmr_data")) {
          dat_agd_contrast <- .unnest_integration(x$network$agd_contrast) %>%
            dplyr::mutate(.sample_size = .data$.sample_size / x$network$n_int)
        } else {
          dat_agd_contrast <- x$network$agd_contrast
        }

        # Only take necessary columns
        dat_agd_contrast <- get_model_data_columns(dat_agd_contrast, regression = x$regression, label = "AgD (contrast-based)")
      } else {
        dat_agd_contrast <- tibble::tibble()
      }

      if (has_ipd(x$network)) {
        dat_ipd <- x$network$ipd
        dat_ipd$.sample_size <- 1

        # Only take necessary columns
        dat_ipd <- get_model_data_columns(dat_ipd, regression = x$regression, label = "IPD")
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

      # Check all variables are present
      regdat <- get_model_data_columns(dat_studies, regression = x$regression, label = "`newdata`")
    }


    # Get number of treatments, number of studies
    ntrt <- nlevels(x$network$treatments)
    nstudy <- nrow(dat_studies)

    # Expand rows for every treatment
    all_trts <- tidyr::expand_grid(.study = dat_studies$.study, .trt = x$network$treatments[-1])
    if (rlang::has_name(dat_studies, ".trt")) dat_studies <- dplyr::select(dat_studies, -".trt")
    dat_studies <- dplyr::left_join(all_trts, dat_studies, by = ".study")

    # Add in .trtclass if defined in network
    if (!is.null(x$network$classes)) {
      dat_studies$.trtclass <- x$network$classes[as.numeric(dat_studies$.trt)]
    }

    # Get model formula and design matrix
    nma_formula <- make_nma_formula(x$regression,
                                    consistency = x$consistency,
                                    classes = !is.null(x$network$classes),
                                    class_interactions = x$class_interactions)

    X_list <- make_nma_model_matrix(nma_formula,
                                    dat_agd_arm = dat_studies,
                                    xbar = x$xbar,
                                    consistency = x$consistency,
                                    classes = !is.null(x$network$classes),
                                    newdata = TRUE)
    X_all <- X_list$X_agd_arm

    # Subset design matrix into EM columns and trt columns
    X_d <- X_all[, grepl("^(\\.trt|\\.contr)[^:]+$", colnames(X_all)), drop = FALSE]
    EM_regex <- "(^\\.trt(class)?.+\\:)|(\\:\\.trt(class)?.+$)"
    X_EM <- X_all[, grepl(EM_regex, colnames(X_all)), drop = FALSE]

    # If there are no EMs (regression model had no interactions with .trt) then
    # just return the treatment effects (not study specific)
    if (ncol(X_EM) == 0) {

      if (!predictive_distribution) {
        re_array <- as.array(as.stanfit(x), pars = "d")
      } else {
        # For predictive distribution, use delta_new instead of d
        re_array <- get_delta_new(x)
      }

      # Set treatments vector
      trtb <- x$network$treatments[-1]
      trta <- x$network$treatments[1]

      if (all_contrasts) {
        re_array <- make_all_contrasts(re_array, trt_ref = nrt)

        contrs <- utils::combn(x$network$treatments, 2)
        trtb <- contrs[2, ]
        trta <- contrs[1, ]

      } else if (!is.null(trt_ref) && trt_ref != nrt) {
        d_ref <- re_array[ , , paste0("d[", trt_ref, "]"), drop = FALSE]
        re_array <- sweep(re_array, 1:2, d_ref, FUN = "-")

        # Add in parameter for network ref trt in place of trt_ref
        re_array[ , , paste0("d[", trt_ref, "]")] <- -d_ref
        d_names <- dimnames(re_array)[[3]]
        d_names[d_names == paste0("d[", trt_ref, "]")] <- paste0("d[", nrt, "]")
        dimnames(re_array)[[3]] <- d_names

        # Reorder parameters
        d_names <- c(paste0("d[", nrt, "]"), d_names[d_names != paste0("d[", nrt, "]")])
        re_array <- re_array[ , , d_names, drop = FALSE]

        trt_reorder <- sort(forcats::fct_relevel(x$network$treatments, trt_ref))
        trtb <- trt_reorder[-1]
        trta <- trt_reorder[1]
      }

      # Fix up parameter names for predictive_distribution = TRUE
      if (predictive_distribution) {
        dimnames(re_array)[[3]] <- gsub("^d\\[", "delta_new[", dimnames(re_array)[[3]])
      }

      if (summary) {
        re_summary <- summary_mcmc_array(re_array, probs = probs)  %>%
          tibble::add_column(.trtb = trtb, .trta = trta, .before = 1)

        out <- list(summary = re_summary, sims = re_array)
      } else {
        out <- list(sims = re_array)
      }
    } else {

      # Which covariates are EMs
      EM_col_names <- stringr::str_remove(colnames(X_EM), EM_regex)
      EM_vars <- get_EM_vars(nma_formula)

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
      rownames(X_EM_d) <- paste0("d[", dat_studies$.trt , "]")

      d_array <- as.array(x, pars = colnames(X_EM_d))
      if (predictive_distribution) {
        # For predictive distribution, use delta_new instead of d
        d_array[ , , colnames(X_d)] <- get_delta_new(x)
      }

      # Linear combination with posterior MCMC array
      re_array <- tcrossprod_mcmc_array(d_array, X_EM_d)

      # Set treatments vector
      trtb <- x$network$treatments[-1]
      trta <- rep(x$network$treatments[1], times = ntrt - 1)

      # Produce all contrasts, if required
      if (all_contrasts) {
        n_contr <- ntrt * (ntrt - 1) / 2
        re_array_all <- array(dim = c(dim(re_array)[1:2],
                                      dplyr::n_distinct(dat_studies$.study) * n_contr),
                              dimnames = purrr::list_modify(dimnames(re_array),
                                                            parameters = paste(rep(1:dplyr::n_distinct(dat_studies$.study), each = n_contr),
                                                                               rep(1:n_contr, times = dplyr::n_distinct(dat_studies$.study)))))
        for (j in unique(dat_studies$.study)) {
          j_select <- dat_studies$.study == j
          j_select_all <- rep(unique(dat_studies$.study) == j, each = n_contr)
          j_contrs_all <- make_all_contrasts(re_array[ , , j_select, drop = FALSE], trt_ref = nrt)
          re_array_all[ , , j_select_all] <- j_contrs_all
          dimnames(re_array_all)[[3]][j_select_all] <- dimnames(j_contrs_all)[[3]]
        }
        re_array <- re_array_all

        contrs <- utils::combn(x$network$treatments, 2)
        trtb <- contrs[2, ]
        trta <- contrs[1, ]
      }

      # Add in study names to parameters
      parnames <- stringr::str_extract(dimnames(re_array)[[3]], "(?<=^d\\[)(.+)(?=\\]$)")
      if (all_contrasts) {
        study_parnames <- paste0("d[", rep(unique(dat_studies$.study), each = n_contr), ": ", parnames, "]")
      } else {
        study_parnames <- paste0("d[", dat_studies$.study, ": ", parnames, "]")
      }
      dimnames(re_array)[[3]] <- study_parnames

      # Rebase against ref_trt if given
      if (!is.null(trt_ref) && trt_ref != nrt) {
        for (j in unique(dat_studies$.study)) {
          j_pars <- dat_studies$.study == j
          d_ref <- re_array[ , , paste0("d[", j, ": ", trt_ref, "]"), drop = FALSE]
          re_array[ , , j_pars] <- sweep(re_array[ , , j_pars, drop = FALSE], 1:2, d_ref, FUN = "-")

          # Add in parameter for network ref trt in place of trt_ref
          re_array[ , , paste0("d[", j, ": ", trt_ref, "]")] <- -d_ref
          d_names <- dimnames(re_array)[[3]]
          d_names[d_names == paste0("d[", j, ": ", trt_ref, "]")] <- paste0("d[", j, ": ", nrt, "]")
          dimnames(re_array)[[3]] <- d_names

          # Reorder parameters
          d_names[j_pars] <- c(paste0("d[", j, ": ", nrt, "]"), d_names[j_pars & d_names != paste0("d[", j, ": ", nrt, "]")])
          re_array <- re_array[ , , d_names, drop = FALSE]

          trt_reorder <- sort(forcats::fct_relevel(x$network$treatments, trt_ref))
          trtb <- trt_reorder[-1]
          trta <- rep(trt_reorder[1], times = ntrt - 1)
        }
      }

      # Fix up parameter names for predictive_distribution = TRUE
      if (predictive_distribution) {
        dimnames(re_array)[[3]] <- gsub("^d\\[", "delta_new[", dimnames(re_array)[[3]])
      }

      # Create summary stats
      if (summary) {
        re_summary <- summary_mcmc_array(re_array, probs = probs) %>%
          tibble::add_column(.study =
                               if (all_contrasts) {
                                 rep(unique(dat_studies$.study), each = n_contr)
                               } else {
                                 dat_studies$.study
                               },
                             .trtb = rep(trtb, times = nstudy),
                             .trta = rep(trta, times = nstudy),
                             .before = 1)
      }

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

      if (summary) {
        out <- list(summary = re_summary, sims = re_array, studies = study_EMs)
      } else {
        out <- list(sims = re_array, studies = study_EMs)
      }
    }
  }

  # Return nma_summary object
  if (summary) {
    class(out) <- "nma_summary"
    attr(out, "xlab") <- if (all_contrasts) "Contrast" else "Treatment"
    attr(out, "ylab") <- get_scale_name(likelihood = x$likelihood,
                                        link = x$link,
                                        measure = "relative",
                                        type = "link")
  }
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

  if (inherits(x, "Matrix")) {
    for (i in 1:nchains) {
      out[ , i, ] <- as.matrix(Matrix::tcrossprod(a[ , i, ], x))
    }
  } else {
    for (i in 1:nchains) {
      out[ , i, ] <- tcrossprod(a[ , i, ], x)
    }
  }

  dimnames(out) <- dimnames_out

  return(out)
}

#' Sample predictive distribution of deltas for new study
#'
#' @param x Fitted RE stan_nma object
#' @param ...
#'
#' @return A 3D MCMC array
#' @noRd
get_delta_new <- function(x, ...) {
  ntrt <- nlevels(x$network$treatments)

  RE_cor <- matrix(0.5, nrow = ntrt - 1, ncol = ntrt - 1)
  diag(RE_cor) <- 1

  draws <- as.matrix(x, pars = c("d", "tau"))
  # Need to make d[] index numeric for rstan::gqs
  colnames(draws)[1:(ntrt - 1)] <- paste0("d[", 1:(ntrt - 1), "]")

  if (ntrt == 2) {
    # Work around rstan::gqs() bug when there is only d[1]
    draws <- draws[, c("d[1]", "d[1]", "tau")]
    colnames(draws) <- c("d[1]", "d[2]", "tau")
    delta_new <- rstan::gqs(stanmodels$predict_delta_new,
                            data = list(nt = 3,
                                        RE_cor = diag(2)),
                            draws = draws,
                            seed = rstan::get_seed(x$stanfit))

    delta_new <- as.array(delta_new)[ , , 1, drop = FALSE]
  } else {
    delta_new <- rstan::gqs(stanmodels$predict_delta_new,
                            data = list(nt = ntrt,
                                        RE_cor = RE_cor),
                            draws = draws,
                            seed = rstan::get_seed(x$stanfit))

    delta_new <- as.array(delta_new)
  }

  # Note: name these d instead of delta_new here so that they can be drop-in
  # replacements in relative_effects(), then fix up later
  dimnames(delta_new)[[3]] <- paste0("d[", levels(x$network$treatments)[-1], "]")
  class(delta_new) <- c("mcmc_array", class(delta_new))

  return(delta_new)
}
