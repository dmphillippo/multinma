#' Compare population covariate distributions
#'
#' Covariate distributions for studies in the network may be compared by overlap
#' effective sample size (via propensity scores) or by (standardised) Euclidean
#' distance between the covariate means. The overlap effective sample size
#' accounts for the full multivariate joint distribution of the covariates,
#' whereas the Euclidean distance between the covariate means only accounts for
#' the location and not the range of variation or joint structure.
#'
#' @param network An `nma_data` network object.
#' @param covariates Character vector of covariate names to compare on.
#' @param method Method to compare distributions, either `"propensity"` to
#' calculate overlap effective sample size with propensity scores, or
#' `"euclidean"` to calculate standardised Euclidean distance between means.
#' @param na_action Optional function (or function name) for handling missing
#'   IPD covariate values, passed to
#'   \code{\link[stats:model.matrix]{model.matrix()}}. The default is to remove
#'   IPD observations with NAs with a warning. Missing values in the AgD will
#'   always give an error.
#'
#' @return A `pop_diff` object, containing a `summary` data frame and
#'   `comparison_matrix` matrix of pairwise comparisons.
#'
#' @importFrom stats glm model.matrix predict weighted.mean sd var weights
#' @export
compare_populations <- function(network,
                                covariates = NULL,
                                method = c("propensity", "euclidean"),
                                na_action) {

  method <- rlang::arg_match(method)

  if (missing(na_action)) na_action <- na_omit_warn

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  # Check for integration call if AgD present
  if (method == "propensity") {
    if ((has_agd_arm(network) || has_agd_contrast(network))) {
      if (is.null(network$int_call)) {
        abort(c('Integration points must be present for method = "propensity".',
                'Set up integration points using `add_integration()` to define the covariate distributions, or set method = "euclidean" to compare means.'))
      } else {
        int_covariates <- names(network$int_call)
      }
    } else {
      int_covariates <- character()
    }
  }

  # Check covariates argument
  if (is.null(covariates)) {
    if (method == "euclidean" || length(int_covariates) < 1) {
      abort('Please provide `covariates` to compare on when method = "euclidean"')
    } else {
      covariates <- int_covariates
      inform(paste0("Comparing on all covariates with integration points: ", paste(covariates, collapse = ", ")))
    }
  } else {
    missing_covs <- setdiff(covariates, int_covariates)
    if (length(missing_covs) > 0) {
      abort(c(paste0("Cannot compare requested covariates missing integration points: ",
                     paste(missing_covs, collapse = ", "),
                   "."), "Set up integration points using `add_integration()`."))
    }
  }

  # ==========================================
  # METHOD: PROPENSITY
  # ==========================================
  if (method == "propensity") {

    # Prepare IPD Data
    if (!isTRUE(nrow(network$ipd) > 0)) {
      abort("IPD must be present when wanting to compare populations using method = `propensity`")
    }

    ipd_df <- as.data.frame(
      lapply(network$ipd[c(covariates, ".study")], function(x) if (is.logical(x)) as.numeric(x) else x)
    )
    ipd_covariate_data <- split(ipd_df[covariates], ipd_df$.study, drop = TRUE)

    # Prepare AGD Data: unnest existing integration points, combine arms within each study
    agd_covariate_data <- if (isTRUE(nrow(network$agd_arm) > 0)) {
      agd_study_labels <- as.character(network$agd_arm$.study)
      arm_data <- lapply(seq_len(nrow(network$agd_arm)), function(i) {
        unnested <- unnest_integration(network$agd_arm[i, , drop = FALSE])
        unnested[covariates]
      })
      names(arm_data) <- agd_study_labels
      # Combine arms that belong to the same study
      lapply(split(arm_data, agd_study_labels), function(arms) do.call(rbind, arms))
    } else {
      list()
    }

    all_data <- c(ipd_covariate_data, agd_covariate_data)
    study_names <- names(all_data)
    n_studies <- length(all_data)

    # Single pass: fit model, compute propensity scores, weights, and ESS per pair
    propensity_scores_list <- list()
    ess_rows <- vector("list", n_studies * (n_studies - 1L) / 2L)
    k <- 0L

    for (i in 1:(n_studies - 1)) {
      for (j in (i + 1):n_studies) {
        s1 <- study_names[i]
        s2 <- study_names[j]

        combined_df <- rbind(
          cbind(all_data[[i]], study_indicator = 1L),
          cbind(all_data[[j]], study_indicator = 0L)
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

        pair_name <- paste(s1, s2, sep = "_vs_")
        propensity_scores_list[[pair_name]] <- combined_df

        k <- k + 1L
        ess_rows[[k]] <- data.frame(
          comparison = pair_name,
          study1 = s1,
          study2 = s2,
          original_n = nrow(combined_df),
          ess = ess,
          ess_percent_of_original = ess / nrow(combined_df) * 100,
          stringsAsFactors = FALSE
        )
      }
    }

    ess_summary <- do.call(rbind, ess_rows)
    ess_summary <- ess_summary[order(ess_summary$ess_percent_of_original, decreasing = TRUE), ]

    # Build symmetric ESS matrix directly
    sorted_names <- sort(study_names)
    n <- length(sorted_names)
    sorted_matrix <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(sorted_names, sorted_names))
    for (r in seq_len(nrow(ess_summary))) {
      s1 <- ess_summary$study1[r]
      s2 <- ess_summary$study2[r]
      val <- ess_summary$ess_percent_of_original[r]
      sorted_matrix[s1, s2] <- val
      sorted_matrix[s2, s1] <- val
    }

    # Detect subnetworks
    g <- igraph::as.igraph(network)
    components <- igraph::components(g)
    treatment_components <- data.frame(
      .trt = names(components$membership),
      subnetwork = components$membership
    )

    study_trt_lookup <- list(network$ipd, network$agd_contrast, network$agd_arm) %>%
      purrr::compact() %>%
      purrr::map_dfr(~ {
        if (all(c(".study", ".trt") %in% colnames(.x))) {
          dplyr::tibble(.study = as.character(.x$.study), .trt = as.character(.x$.trt))
        } else {
          NULL
        }
      }) %>%
      dplyr::distinct()

    study_components <- study_trt_lookup %>%
      dplyr::left_join(treatment_components, by = ".trt") %>%
      dplyr::select(-".trt") %>%
      dplyr::distinct(.data$.study, .data$subnetwork)

    output_list <- list(
      propensity_scores = propensity_scores_list,
      summary = ess_summary,
      full_matrix = sorted_matrix
    )

    if (max(study_components$subnetwork) == 2) {
      sub1 <- dplyr::filter(study_components, .data$subnetwork == 1)$.study
      sub2 <- dplyr::filter(study_components, .data$subnetwork == 2)$.study
      output_list$subnetwork_matrix <- sorted_matrix[
        rownames(sorted_matrix) %in% sub1,
        colnames(sorted_matrix) %in% sub2,
        drop = FALSE
      ]
    }

    return(output_list)

  }

  # ==========================================
  # METHOD: EUCLIDEAN
  # ==========================================
  if (method == "euclidean") {

    # Validation
    if (isTRUE(nrow(network$agd_contrast) > 0)) {
      # Check if the specific column ".sample_size" is MISSING
      if (!".sample_size" %in% colnames(network$agd_contrast)) {
        abort("Aggregate contrast data must contain a '.sample_size' column.")
      }
    }

    # Process IPD: Convert logical to numeric and average by study
    if (nrow(network$ipd) > 0) {
      ipd_covariate_data <- network$ipd
      ipd_covariate_data[covariates] <- lapply(ipd_covariate_data[covariates], function(x) {
        if (is.logical(x)) as.numeric(x) else x
      })
      ipd_summary <- ipd_covariate_data %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::summarise(
          total_n = dplyr::n(),
          dplyr::across(dplyr::all_of(covariates), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
          .groups = "drop"
        )
    }

    # Helper Function: Extract AGD Means
    extract_agd_means <- function(agd_df) {
      if (nrow(agd_df) == 0) return(NULL)

      df <- agd_df[, c(".study", ".sample_size"), drop = FALSE]
      retained_covariates <- c()

      for (cov in covariates) {
        distr_obj <- network$int_call[[cov]]
        args <- distr_obj

        is_binary_dist <- "prob" %in% names(args)

        add_covariate <- TRUE

        if (is_binary_dist) {
          # BINARY LOGIC (Bernoulli/Binomial)
          prob_col <- as.character(args$prob)

          if (prob_col %in% colnames(agd_df)) {
            p <- agd_df[[prob_col]]
            df[[paste0(cov, "_mean")]] <- p
            # Auto-calculate SD for binary: sqrt(p * (1-p))
            df[[paste0(cov, "_sd")]] <- sqrt(p * (1 - p))
          } else {
            warn(glue::glue("Probability column '{prob_col}' for covariate '{cov}' not found in AgD. Dropped."))
            add_covariate <- FALSE
          }

        } else {

          mean_col <- as.character(args$mean)
          sd_col <- as.character(args$sd)

        # Handle Means
        if (mean_col %in% colnames(agd_df)) {
          df[[paste0(cov, "_mean")]] <- agd_df[[mean_col]]
        } else if (cov %in% colnames(agd_df)) {
          df[[paste0(cov, "_mean")]] <- agd_df[[cov]]
        } else {
          warn(glue::glue("Mean for covariate '{cov}' not found in AgD. Covariate dropped."))
          add_covariate <- FALSE
        }

        # Handle SDs
        if (add_covariate) {
          if (sd_col %in% colnames(agd_df)) {
            # Explicit SD column exists -> Use it
            df[[paste0(cov, "_sd")]] <- agd_df[[sd_col]]

          } else {
              # Continuous variable missing SD (or percentage > 1) -> Drop it
              warn(glue::glue("SD column '{sd_col}' for covariate '{cov}' not found in AgD. Dropped (Continuous variable requires explicit SD)."))
              df[[paste0(cov, "_mean")]] <- NULL
              add_covariate <- FALSE
            }
          }
        }

        # --- Final Decision ---
        if (add_covariate) {
          retained_covariates <- c(retained_covariates, cov)
        } else {
          # Cleanup if we failed halfway through
          df[[paste0(cov, "_mean")]] <- NULL
          df[[paste0(cov, "_sd")]] <- NULL
        }
      }
      # Update the 'covariates' list to exclude dropped ones
      assign("covariates", retained_covariates, envir = parent.env(environment()))
      return(df)
    }

    # Extract AGD Data
    agd_contrast_means <- extract_agd_means(network$agd_contrast)
    agd_arm_means <- extract_agd_means(network$agd_arm)
    agd_all <- dplyr::bind_rows(agd_contrast_means, agd_arm_means)

    # Check for NAs in extracted AGD
    idx <- which(is.na(agd_all), arr.ind = TRUE)
    if (nrow(idx)) {
      rows <- idx[, "row"]
      cols <- idx[, "col"]
      studies <- if (".study" %in% names(agd_all)) agd_all$.study[rows] else rownames(agd_all)[rows]
      vars <- colnames(agd_all)[cols]
      miss <- unique(data.frame(study = studies, variable = vars, stringsAsFactors = FALSE))
      lines_by_var <- tapply(miss$study, miss$variable, function(s) paste(unique(s), collapse = ", "))
      abort(paste0(
        "AgD covariate inputs contain missing values:\n",
        paste(" \u2022 ", names(lines_by_var), " missing in studies: ", unname(lines_by_var), collapse = "\n"),
        "\nPlease remove these variables from `covariates`"
      ))
    }

    # Weighted Summary of AGD
    agd_summary <- agd_all %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise(
        total_n = sum(.data$.sample_size, na.rm = TRUE),
        !!!setNames(unlist(lapply(covariates, function(cov) {
          list(
            rlang::expr(weighted.mean(!!rlang::sym(paste0(cov, "_mean")), w = .data$.sample_size, na.rm = TRUE)),
            rlang::expr(sqrt(weighted.mean((!!rlang::sym(paste0(cov, "_sd")))^2, w = .data$.sample_size, na.rm = TRUE)))
          )
        }), recursive = FALSE), unlist(lapply(covariates, function(cov) {
          c(paste0(cov, "_mean"), paste0(cov, "_sd"))
        }))),
        .groups = "drop"
      )

    ipd_summary$source <- "IPD"
    agd_summary$source <- "AGD"
    all_summary <- dplyr::bind_rows(ipd_summary, agd_summary)

    # Identify Subnetworks
    g <- igraph::as.igraph(network)
    components <- igraph::components(g)
    treatment_components <- data.frame(
      .trt = names(components$membership),
      subnetwork = components$membership
    )
    study_trt_lookup <- list(network$ipd, network$agd_contrast, network$agd_arm) %>%
      purrr::compact() %>%
      purrr::map_dfr(~ {
        cols <- colnames(.x)
        if (all(c(".study", ".trt") %in% cols)) {
          dplyr::tibble(.study = as.character(.x$.study), .trt = as.character(.x$.trt))
        } else {
          NULL
        }
      }) %>%
      dplyr::distinct()

    study_components <- study_trt_lookup %>%
      dplyr::left_join(treatment_components, by = ".trt") %>%
      dplyr::select(-".trt") %>%
      dplyr::distinct(.data$.study, .data$subnetwork)

    all_summary <- dplyr::left_join(all_summary, study_components, by = ".study")

    sub1 <- dplyr::filter(all_summary, .data$subnetwork == 1)
    sub2 <- dplyr::filter(all_summary, .data$subnetwork == 2)

    dist_matrix <- matrix(NA, nrow = nrow(sub1), ncol = nrow(sub2), dimnames = list(sub1$.study, sub2$.study))
    dist_matrix_full <- matrix(NA, nrow = nrow(all_summary), ncol = nrow(all_summary), dimnames = list(all_summary$.study, all_summary$.study))

    # Calculate Distances (Subnetwork 1 vs 2)
    if (nrow(sub2) < 1) {
      for (i in seq_len(nrow(sub1))) {
        for (j in seq_len(nrow(sub2))) {
          vec1 <- as.numeric(sub1[i, paste0(covariates, "_mean")])
          vec2 <- as.numeric(sub2[j, paste0(covariates, "_mean")])
          sd1 <- as.numeric(sub1[i, paste0(covariates, "_sd")])
          sd2 <- as.numeric(sub2[j, paste0(covariates, "_sd")])
          pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
          valid <- !is.na(vec1) & !is.na(vec2) & !is.na(pooled_sd) & pooled_sd > 0
          diff_scaled <- (vec1[valid] - vec2[valid]) / pooled_sd[valid]
          dist_matrix[i, j] <- sqrt(sum(diff_scaled^2))
        }
      }
    }

    # Calculate Distances (Full Matrix)
    for (i in seq_len(nrow(all_summary))) {
      for (j in seq_len(nrow(all_summary))) {
        if (i == j) next
        vec1 <- as.numeric(all_summary[i, paste0(covariates, "_mean")])
        vec2 <- as.numeric(all_summary[j, paste0(covariates, "_mean")])
        sd1 <- as.numeric(all_summary[i, paste0(covariates, "_sd")])
        sd2 <- as.numeric(all_summary[j, paste0(covariates, "_sd")])
        pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
        valid <- !is.na(vec1) & !is.na(vec2) & !is.na(pooled_sd) & pooled_sd > 0
        diff_scaled <- (vec1[valid] - vec2[valid]) / pooled_sd[valid]
        dist_matrix_full[i, j] <- sqrt(sum(diff_scaled^2))
      }
    }

    # Return Euclidean Results
    if (nrow(sub2) < 1) {
      return(list(summary = all_summary, distance_matrix = dist_matrix_full))
    } else {
      return(list(summary = all_summary, distance_matrix = dist_matrix, distance_matrix_full = dist_matrix_full))
    }
  }
}

na_omit_warn <- function(object, ...) {
  out <- stats::na.omit(object, ...)
  if (length(out) < length(object)) warn("Observations with missing values removed.")
  return(out)
}

