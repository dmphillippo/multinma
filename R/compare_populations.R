#' Compare population covariate distributions
#'
#' Covariate distributions for studies in the network may be compared by overlap
#' effective sample size (via propensity scores) or by (standardised) Euclidean
#' distance between the covariate means. The overlap effective sample size
#' accounts for the full multivariate joint distribution of the covariates,
#' whereas the Euclidean distance between the covariate means only accounts for
#' the location and not the range of variation or joint structure.
#'
#' # Aggregate data setup
#' When aggregate data are present in the network, some setup is necessary to
#' use this function. For `method = "propensity"`, the network must have
#' integration points present; specify these using [add_integration()]. For
#' `method = "euclidean"`, the network must either have integration points
#' present, or the mean and standard deviation of covariates can be provided in
#' the original input data when setting up the network. In the latter case, the
#' covariate means are assumed to be provided in the named `<covariate>` column,
#' and the corresponding standard deviation in a column named `<covariate>_sd`.
#' For example, if the covariate is `age`, mean age should be provided in the
#' `age` column, and the standard devation of age in `age_sd` column.
#'
#' @param network An `nma_data` network object.
#' @param covariates Character vector of covariate names to compare on.
#' @param method Method to compare distributions, either `"propensity"` to
#' calculate overlap effective sample size with propensity scores, or
#' `"euclidean"` to calculate standardised Euclidean distance between means.
#'
#' @return A `pop_comp` object, containing a `summary` data frame and
#'   `comparison_matrix` matrix of pairwise comparisons.
#'   When `method = "propensity"`, a list `propensity_scores` of data frames of
#'   fitted propensity scores will also be included.
#'
#' @importFrom stats glm model.matrix predict weighted.mean sd var weights
#' @export
compare_populations <- function(network,
                                covariates = NULL,
                                method = c("propensity", "euclidean")) {

  method <- rlang::arg_match(method)

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  # Check for integration call if AgD present
  if ((has_agd_arm(network) || has_agd_contrast(network))) {
    if (is.null(network$int_call)) {
      if (method == "propensity") {
        abort(c('Integration points must be present for method = "propensity".',
                'Set up integration points using `add_integration()` to define the covariate distributions, or set method = "euclidean" to compare means.'))
      } else {
        int_covariates <- character()
      }
    } else {
      int_covariates <- names(network$int_call)
    }
  } else {
    int_covariates <- character()
  }

  # Need AgD sample sizes
  if (!has_agd_sample_size(network)) {
    abort(c("AgD sample sizes must be available.",
            "Provide these using the `sample_size` argument in set_agd_*()."))
  }

  # Check covariates argument
  if (is.null(covariates)) {
    if (method == "euclidean" && length(int_covariates) < 1) {
      abort('Provide `covariates` to compare on when method = "euclidean"')
    } else if (method == "propensity" && length(int_covariates) < 1) {
      abort('Provide `covariates` to compare on')
    } else {
      covariates <- int_covariates
      inform(paste0("Comparing on all covariates with integration points: ", paste(covariates, collapse = ", ")))
    }
  } else if ((has_agd_arm(network) || has_agd_contrast(network)) && !is.null(network$int_call)) {
    missing_covs <- setdiff(covariates, int_covariates)
    if (length(missing_covs) > 0) {
      abort(c(paste0("Cannot compare requested covariates missing integration points: ",
                     paste(missing_covs, collapse = ", "),
                     "."), "Set up integration points using `add_integration()`."))
    }
  }

  # Get covariates
  cov_formula <- as.formula(paste0("~", paste(covariates, collapse = " + ")))
  if (has_ipd(network)) {
    dat_ipd <- network$ipd

    if (!all(covariates %in% names(dat_ipd))) {
      abort(paste0("Covariates not found in IPD:", paste(setdiff(covariates, names(dat_ipd)), sep = ", ")))
    }

    complete <- complete.cases(dat_ipd[, covariates])
    if (any(!complete)) {
      nmiss <- sum(!complete)
      warn(glue::glue("Removed {nmiss} observation{if (nmiss > 1) 's' else ''} with missing covariate values from IPD."))
      dat_ipd <- dat_ipd[complete, ]
    }

    withCallingHandlers(
      ipd_covs <- as.data.frame(model.matrix(cov_formula, dat_ipd)[, -1]),
      error = ~abort(paste0("Failed to get IPD covariate data.\n", .)))
    ipd_study <- dat_ipd$.study
  } else {
    ipd_covs <- ipd_study <- NULL
  }

  if (has_agd_arm(network)) {
    if (!is.null(network$int_call)) {
      # Resample integration points with n_int = .sample_size
      ds <- purrr::map(rlang::list2(!!! network$int_call), rlang::eval_tidy)

      dat_agd_arm <- network$agd_arm %>%
        dplyr::select(-dplyr::starts_with(".int_")) %>%
        dplyr::group_by(.data$.study, .data$.trt) %>%
        dplyr::group_modify(~rlang::exec(add_integration, x = .x,
                                         !!! ds,
                                         n_int = .x$.sample_size,
                                         cor = network$int_cor)) %>%
        .unnest_integration()
    } else {
      dat_agd_arm <- network$agd_arm
    }

    if (!all(covariates %in% names(dat_agd_arm))) {
      abort(paste0("Covariates not found in AgD (arm-based): ", paste(setdiff(covariates, names(dat_agd_arm)), collapse = ", ")))
    }

    complete <- complete.cases(dat_agd_arm[, covariates])
    if (any(!complete)) {
      nmiss <- nrow(dplyr::distinct(dat_agd_arm[!complete, ], .data$.study, .data$.trt))
      warn(glue::glue("Removed {nmiss} study arm{if (nmiss > 1) 's' else ''} with missing covariate values from AgD (arm-based)."))

      dat_agd_arm <- dat_agd_arm[complete, ]
    }

    withCallingHandlers(
      agd_arm_covs <- as.data.frame(model.matrix(cov_formula, dat_agd_arm)[, -1]),
      error = ~abort(paste0("Failed to get Agd (arm-based) covariate data.\n", .)))
    agd_arm_study <- dat_agd_arm$.study
  } else {
    agd_arm_covs <- agd_arm_study <- NULL
  }

  if (has_agd_contrast(network)) {
    if (!is.null(network$int_call)) {
      # Resample integration points with n_int = .sample_size
      ds <- purrr::map(rlang::list2(!!! network$int_call), rlang::eval_tidy)

      dat_agd_contrast <- network$agd_contrast %>%
        dplyr::select(-dplyr::starts_with(".int_")) %>%
        dplyr::group_by(.data$.study, .data$.trt) %>%
        dplyr::group_modify(~rlang::exec(add_integration, x = .x,
                                         !!! ds,
                                         n_int = .x$.sample_size,
                                         cor = network$int_cor)) %>%
        .unnest_integration()

    } else {
      dat_agd_contrast <- network$agd_contrast
    }

    if (!all(covariates %in% names(dat_agd_contrast))) {
      abort(paste0("Covariates not found in AgD (contrast-based):", paste(setdiff(covariates, names(dat_agd_contrast)), sep = ", ")))
    }

    complete <- complete.cases(dat_agd_contrast[, covariates])
    if (any(!complete)) {
      nmiss <- nrow(dplyr::distinct(dat_agd_contrast[!complete, ], .data$.study, .data$.trt))
      warn(glue::glue("Removed {nmiss} study arm{if (nmiss > 1) 's' else ''} with missing covariate values from AgD (contrast-based)."))

      dat_agd_contrast <- dat_agd_contrast[complete, ]
    }

    withCallingHandlers(
      agd_contrast_covs <- as.data.frame(model.matrix(cov_formula, dat_agd_contrast)[, -1]),
      error = ~abort(paste0("Failed to get Agd (contrast-based) covariate data.\n", .)))
    agd_contrast_study <- dat_agd_contrast$.study
  } else {
    agd_contrast_covs <- agd_contrast_study <- NULL
  }


  if (method == "propensity") {

    dat_list <- c(if (has_ipd(network)) split(ipd_covs, ipd_study, drop = TRUE) else NULL,
                  if (has_agd_arm(network)) split(agd_arm_covs, agd_arm_study, drop = TRUE) else NULL,
                  if (has_agd_contrast(network)) split(agd_contrast_covs, agd_contrast_study, drop = TRUE) else NULL)

    studies <- network$studies
    n_studies <- length(dat_list)

    # Single pass: fit model, compute propensity scores, weights, and ESS per pair
    propensity_scores <- list()
    ess_rows <- vector("list", n_studies * (n_studies - 1L) / 2L)
    k <- 1L

    for (i in 1:(n_studies - 1)) {
      for (j in (i + 1):n_studies) {
        s1 <- studies[i]
        s2 <- studies[j]

        combined_df <- rbind(
          cbind(dat_list[[s1]], study_indicator = 1L),
          cbind(dat_list[[s2]], study_indicator = 0L)
        )

        model <- stats::glm(study_indicator ~ ., data = combined_df, family = "binomial")
        ps <- stats::predict(model, newdata = combined_df, type = "response")

        combined_df$propensity_score <- ps
        combined_df$overlap_weight <- ifelse(
          combined_df$study_indicator == 1L,
          1 / ps,
          1 / (1 - ps)
        )

        w <- combined_df$overlap_weight
        ess <- sum(w)^2 / sum(w^2)

        pair_name <- paste(s1, s2, sep = " vs. ")
        propensity_scores[[pair_name]] <- combined_df

        ess_rows[[k]] <- dplyr::tibble(
          comparison = pair_name,
          study1 = s1,
          study2 = s2,
          original_n = nrow(combined_df),
          ess = ess,
          ess_percent = ess / nrow(combined_df) * 100
        )
        k <- k + 1L
      }
    }

    ess_summary <- dplyr::bind_rows(ess_rows)
    ess_summary <- ess_summary[order(ess_summary$ess_percent, decreasing = TRUE), ]

    # Build symmetric ESS matrix directly
    sorted_matrix <- matrix(100, nrow = n_studies, ncol = n_studies, dimnames = list(studies, studies))
    for (r in seq_len(nrow(ess_summary))) {
      s1 <- ess_summary$study1[r]
      s2 <- ess_summary$study2[r]
      val <- ess_summary$ess_percent[r]
      sorted_matrix[s1, s2] <- val
      sorted_matrix[s2, s1] <- val
    }

    out <- list(
      summary = ess_summary,
      comparison_matrix = sorted_matrix,
      propensity_scores = propensity_scores
    )

  } else if (method == "euclidean") {

    cov_sd <- paste0(covariates, "_sd")

    if (has_ipd(network)) {
      ipd_summary <- dplyr::mutate(ipd_covs, .study = ipd_study) %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::summarise(
          sample_size = dplyr::n(),
          dplyr::across(dplyr::all_of(covariates), list(mean = mean, sd = sd)),
          .groups = "drop"
        )
    } else {
      ipd_summary <- NULL
    }

    if (has_agd_arm(network) || has_agd_contrast(network)) {
      dat_agd_all <- dplyr::bind_rows(network$agd_arm, network$agd_contrast)
      agd_covs_all <- dplyr::bind_rows(agd_arm_covs, agd_contrast_covs)
      agd_study_all <- c(agd_arm_study, agd_contrast_study)

      if (!is.null(network$int_call)) {
        # With integration points available, use these to get summaries - will be identical to using reported stats
        agd_summary <- dplyr::mutate(agd_covs_all, .study = agd_study_all) %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::summarise(
            sample_size = dplyr::n(),
            dplyr::across(dplyr::all_of(covariates), list(mean = mean, sd = sd)),
            .groups = "drop"
          )
      } else {
        # One row per study arm, reported summary stats.
        # Assume covariate column provides mean, covariate_sd column available for sd
        cov_sd <- paste0(covariates, "_sd")
        miss_sd <- setdiff(cov_sd, names(dat_agd_all))
        if (length(miss_sd)) {
          abort(glue::glue("Standard deviation columns not found in AgD: ",
                           glue::glue_collapse(miss_sd, sep = ", ", last = " and ")))
        }
        if (any(is.na(dat_agd_all[, cov_sd])) || any(is.infinite(unlist(dat_agd_all[, cov_sd])))) {
          abort("Missing or infinite values for covariate standard deviations.")
        }

        agd_summary <- dplyr::mutate(agd_covs_all, .study = agd_study_all) %>%
          dplyr::bind_cols(dat_agd_all[, c(".sample_size", cov_sd)]) %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::summarise(
            sample_size = sum(.data$.sample_size),
            dplyr::across(dplyr::all_of(covariates), ~weighted.mean(., .data$.sample_size), .names = "{.col}_mean"),
            dplyr::across(dplyr::all_of(cov_sd), ~sqrt(weighted.mean(.^2, .data$.sample_size - 1)), .names = "{.col}")
          )
      }
    } else {
      agd_summary <- NULL
    }

    all_summary <- dplyr::bind_rows(ipd_summary, agd_summary) %>%
      dplyr::arrange(.data$.study)

    # Get overall standardising sd for each covariate
    ssd <- all_summary %>%
      dplyr::ungroup() %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(cov_sd), ~sqrt(weighted.mean(.^2, .data$sample_size - 1)),
                                     .names = "{.col}")) %>%
      unlist()

    dist_matrix <- matrix(0, nrow = nrow(all_summary), ncol = nrow(all_summary),
                          dimnames = list(all_summary$.study, all_summary$.study))

    # Calculate Distances
    cov_mean <- paste0(covariates, "_mean")
    cov_sd <- paste0(covariates, "_sd")
    for (i in seq_len(nrow(all_summary))) {
      for (j in seq_len(nrow(all_summary))) {
        if (i == j) next
        vec1 <- all_summary[i, cov_mean]
        vec2 <- all_summary[j, cov_mean]
        diff_scaled <- (vec1 - vec2) / ssd
        dist_matrix[i, j] <- sqrt(sum(diff_scaled^2))
      }
    }

    # Summary data frame
    summary_df <- dplyr::as_tibble(dist_matrix, rownames = "study1") %>%
      tidyr::pivot_longer(!"study1", names_to = "study2", values_to = "distance") %>%
      dplyr::mutate(comparison = paste(.data$study1, .data$study2, sep = " vs. "),
                    study1 = factor(.data$study1, levels = levels(network$studies)),
                    study2 = factor(.data$study2, levels = levels(network$studies))) %>%
      # Add total sample size
      dplyr::left_join(dplyr::transmute(all_summary, .data$.study, ss1 = .data$sample_size),
                       by = dplyr::join_by(x$study1 == y$.study)) %>%
      dplyr::left_join(dplyr::transmute(all_summary, .data$.study, ss2 = .data$sample_size),
                       by = dplyr::join_by(x$study2 == y$.study)) %>%
      dplyr::mutate(sample_size = ss1 + ss2) %>%
      dplyr::select(-"ss1", -"ss2") %>%
      dplyr::relocate("comparison", "study1", "study2", "sample_size", dplyr::everything()) %>%
      dplyr::filter(which(network$studies == .data$study1) <
                      which(network$studies == .data$study2)) %>%
      dplyr::arrange(distance)


    out <- list(summary = summary_df,
                comparison_matrix = dist_matrix)
  }

  # Detect subnetworks
  g <- igraph::as.igraph(network, collapse = FALSE)
  comps <- igraph::decompose(g)
  components <- purrr::imap_dfr(comps,
    ~dplyr::tibble(.study = unique(igraph::edge_attr(.x, ".study")), component = .y))

  # Common outputs
  out$components <- components
  out$method <- method
  out$covariates <- covariates

  class(out) <- c("pop_comp", class(out))
  return(out)
}


#' @param x A `pop_comp` object produced by `compare_populations()`
#' @param order String, the order in which to display the comparison summaries.
#'   Either `"decreasing"`, to list most similar studies first (the default), or
#'   `"increasing"` to list most dissimilar studies first.
#' @param simplify Logical, should the output be simplified to only show
#'   comparisons between subnetworks (`TRUE`, default), or between every study
#'   in the network (`FALSE`)? If a connected network was provided, then
#'   `simplify` is always `FALSE`.
#' @param ... Additional arguments (unused)
#' @param digits Number of digits to print in the summary, default 2
#' @param n Number of rows to show in the table of study summaries
#'
#' @export
#' @rdname compare_populations
print.pop_comp <- function(x,
                           order = c("decreasing", "increasing"),
                           simplify = TRUE, ..., digits = 2L, n = 10) {

  order <- rlang::arg_match(order)
  if (!rlang::is_bool(simplify))
    abort("`simplify` must be TRUE or FALSE.")
  if (!rlang::is_integerish(digits, n = 1, finite = TRUE) || digits < 0)
    abort("`digits` must be a single non-negative integer.")
  if (!rlang::is_integerish(x = n, n = 1) || n < 1)
    abort("`n` must be a single positive integer.")


  m <- switch(x$method,
              propensity = "propensity score overlap",
              euclidean = "standardised Euclidean distance")
  cglue("Compared populations using {m}, based on the following covariates: ",
        "{glue::glue_collapse(x$covariates, sep = ', ', last = ' and ')}.")
  cat("\n")


  sec_header(glue::glue("Pairwise comparisons in {order} order of similarity"))

  ncomp <- max(x$components$component)
  if (ncomp > 2) {
    cglue(subtle("Subnetwork shown in brackets next to study name."))
  }

  x_sum <- x$summary

  if (order == "decreasing") {
    if (x$method == "propensity") {
      x_sum <- dplyr::arrange(x_sum, dplyr::desc(.data$ess_percent))
    } else if (x$method == "euclidean") {
      x_sum <- dplyr::arrange(x_sum, distance)
    }
  } else {
    if (x$method == "propensity") {
      x_sum <- dplyr::arrange(x_sum, .data$ess_percent)
    } else if (x$method == "euclidean") {
      x_sum <- dplyr::arrange(x_sum, dplyr::desc(distance))
    }
  }

  num_col <- setdiff(names(x_sum)[purrr::map_lgl(x_sum, is.numeric)], c("original_n", "sample_size"))

  if (simplify && ncomp > 1) {
    comp_lookup <- setNames(x$components$component, x$components$.study)

    # Add component numbers if more than 2 components
    if (ncomp > 2) {
      x_sum$comparison <- paste0(x_sum$study1, " (", comp_lookup[x_sum$study1], ") vs. ",
                                 x_sum$study2, " (", comp_lookup[x_sum$study2], ")")
    }

    # Show only comparisons across subnetworks
    x_sum <- dplyr::filter(x_sum, comp_lookup[study1] != comp_lookup[study2])
  }

  x_sum <- dplyr::mutate_at(x_sum, num_col, ~round(., digits)) %>%
    dplyr::select(-"study1", -"study2") %>%
    as.data.frame()

  names(x_sum)[names(x_sum) == "comparison"] <- ""

  print(head(x_sum, n = n), row.names = FALSE)

  if (nrow(x_sum) > n) {
    cglue(subtle(" ... plus {nrow(x_sum) - n} more comparisons"))
  }

  cat("\n")
  sec_header("Matrix of pairwise comparisons")

  mat <- round(x$comparison_matrix, digits)

  if (simplify && ncomp > 1) {
    for (c1 in 1:(ncomp - 1)) for (c2 in 2:ncomp) {
      cat("Subnetwork ", c2, " vs. ", c1, ":\n", sep = "")
      s1 <- which(comp_lookup == c1)
      s2 <- which(comp_lookup == c2)
      print(mat[s1, s2])
      cat("\n")
    }
  } else {
    print(mat)
  }

}

