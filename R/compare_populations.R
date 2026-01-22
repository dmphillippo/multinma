#' Comparing population covariates
#'
#' Runs distance tests between multivariate distributions across studies.
#'
#' @param network An `nma_data` object.
#' @param covariates Character vector of covariate names. Cannot accept categorical variables.
#' @param method Character; "euclidean", or "propensity". Determines the distance calculation method.
#'
#' @return A list with a `summary` dataframe and `distance_matrix`.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select rename filter group_by summarise mutate left_join bind_rows distinct across all_of tibble n
#' @importFrom tidyr separate pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom rlang abort expr sym
#' @importFrom stats glm model.matrix predict weighted.mean sd var weights
#' @importFrom utils type.convert
#' @export
compare_populations <- function(network,
                                covariates = NULL,
                                method = c("euclidean", "propensity")) {

  method <- tryCatch({
    if (missing(method)) stop()
    match.arg(method, choices = valid_choices)
  }, error = function(e) {
    rlang::abort("Argument `method` must be either 'euclidean' or 'propensity'.")
  })

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  # Check for Integration Code
  if (is.null(network$integration_code)) {
    abort("The network object does not contain integration code. Please run `add_integration()` on your network first to define the covariate distributions.")
  }
  avail_covariates <- names(network$integration_code)

  # Checks for covariates argument
  if (is.null(covariates)) {
    covariates <- avail_covariates
    message(paste0("Using covariates defined in integration code: ", paste(covariates, collapse = ", ")))
  } else {
    missing_covs <- setdiff(covariates, avail_covariates)
    if (length(missing_covs) > 0) {
      abort(paste0("The following covariates are NOT defined in the network's integration code: ",
                   paste(missing_covs, collapse = ", "),
                   ". Please define them using `add_integration()` or remove them."))
    }
  }

  # ==========================================
  # METHOD: PROPENSITY
  # ==========================================
  if (method == "propensity") {

    # Prepare IPD Data
    if (isTRUE(nrow(network$ipd) > 0)) {
      ipd_covariate_data <- network$ipd
      columns_to_keep <- c(covariates, ".study")
      ipd_covariate_data <- ipd_covariate_data[columns_to_keep]
      ipd_covariate_data <- as.data.frame(
        lapply(ipd_covariate_data, function(x) {
          if (is.logical(x)) as.numeric(x) else x
        })
      )
      ipd_covariate_data <- split(ipd_covariate_data, ipd_covariate_data$.study, drop = TRUE)
      ipd_covariate_data <- lapply(ipd_covariate_data, function(x) {
        x$.study <- NULL
        return(x)
      })
    } else {
      abort("IPD must be present when wanting to compare populations using method = `propensity`")
    }

    # Prepare AGD Data (via Integration)
    if (isTRUE(nrow(network$agd_arm) > 0)) {
      studies <- network$agd_arm$.study
      arm_indices <- seq_len(nrow(network$agd_arm))

      agd_arm_networks_list <- lapply(arm_indices, function(i) {
        new_network <- network
        new_network$agd_arm <- new_network$agd_arm[i, , drop = FALSE]
        return(new_network)
      })
      names(agd_arm_networks_list) <- studies

      agd_arm_networks_list <- lapply(agd_arm_networks_list, function(net) {
        total_sample_size <- sum(net$agd_arm$.sample_size)
        integration_distr_objects <- lapply(network$integration_code, eval)
        other_args <- list(x = net, n_int = total_sample_size)
        all_args <- c(other_args, integration_distr_objects)

        net_integrated <- do.call(add_integration, all_args)
        unnested_agd <- unnest_integration(net_integrated$agd_arm)
        unnested_agd <- unnested_agd[covariates]
        return(unnested_agd)
      })

      agd_arm_covariate_data <- split(agd_arm_networks_list, names(agd_arm_networks_list))
      agd_arm_covariate_data <- lapply(agd_arm_covariate_data, function(sub_list) {
        do.call(rbind, sub_list)
      })
    }

    # Combine Data
    if (isTRUE(nrow(network$ipd) > 0)) all_data <- ipd_covariate_data
    if (isTRUE(nrow(network$agd_arm) > 0)) all_data <- c(all_data, agd_arm_covariate_data)

    # Run Pairwise Logistic Regressions
    regression_results <- list()

    for (i in 1:(length(all_data) - 1)) {
      for (j in (i + 1):length(all_data)) {
        study1_df <- all_data[[i]]
        study2_df <- all_data[[j]]
        combined_df <- rbind(study1_df, study2_df)
        combined_df$study_indicator <- c(rep(1, nrow(study1_df)), rep(0, nrow(study2_df)))

        model <- glm(study_indicator ~ ., data = combined_df, family = "binomial")
        pair_name <- paste(names(all_data)[i], names(all_data)[j], sep = "_vs_")
        regression_results[[pair_name]] <- model
      }
    }

    # Calculate Propensity Scores & Weights
    propensity_scores_list <- list()
    for (pair_name in names(regression_results)) {
      model <- regression_results[[pair_name]]
      study_names <- strsplit(pair_name, "_vs_")[[1]]

      study1_df <- all_data[[study_names[1]]]
      study1_df$study_indicator <- 1
      study2_df <- all_data[[study_names[2]]]
      study2_df$study_indicator <- 0

      combined_df <- rbind(study1_df, study2_df)
      propensity_scores <- predict(model, newdata = combined_df, type = "response")

      combined_df$propensity_score <- propensity_scores
      combined_df$ate_weight <- ifelse(
        combined_df$study_indicator == 1,
        1 / combined_df$propensity_score,
        1 / (1 - combined_df$propensity_score)
      )
      propensity_scores_list[[pair_name]] <- combined_df
    }

    # Calculate Effective Sample Size (ESS)
    ess_summary <- data.frame(
      comparison = character(),
      original_n = integer(),
      ess = numeric(),
      ess_reduction_percent = numeric(),
      stringsAsFactors = FALSE
    )

    for (pair_name in names(propensity_scores_list)) {
      df <- propensity_scores_list[[pair_name]]
      sum_of_weights <- sum(df$ate_weight)
      sum_of_squared_weights <- sum(df$ate_weight^2)
      effective_sample_size <- (sum_of_weights^2) / sum_of_squared_weights

      original_n <- nrow(df)
      ess_percent <- (effective_sample_size / original_n) * 100

      ess_summary <- rbind(ess_summary, data.frame(
        comparison = pair_name,
        original_n = original_n,
        ess = effective_sample_size,
        ess_percent_of_original = ess_percent
      ))
    }
    ess_summary <- ess_summary[order(ess_summary$ess_percent_of_original, decreasing = TRUE), ]

    # Format Output Matrix (Long to Wide)
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
      dplyr::select(-.trt) %>%
      dplyr::distinct(.study, subnetwork)

    df_long <- ess_summary %>%
      separate(comparison, into = c("item1", "item2"), sep = "_vs_") %>%
      select(item1, item2, value = ess_percent_of_original)

    df_symmetric <- df_long %>%
      rename(item2 = item1, item1 = item2, value = value)

    df_all <- rbind(df_long, df_symmetric)
    final_matrix_df <- df_all %>%
      pivot_wider(names_from = item2, values_from = value)

    final_matrix <- final_matrix_df %>%
      column_to_rownames(var = "item1") %>%
      as.matrix()

    sorted_names <- sort(rownames(final_matrix))
    sorted_matrix <- final_matrix[sorted_names, sorted_names]

    # Filter by subnetwork if applicable
    if (max(study_components$subnetwork) == 2) {
      sub1 <- dplyr::filter(study_components, subnetwork == 1)
      sub2 <- dplyr::filter(study_components, subnetwork == 2)
      rows_to_keep <- rownames(sorted_matrix) %in% sub1$.study
      cols_to_keep <- colnames(sorted_matrix) %in% sub2$.study
      filtered_matrix <- sorted_matrix[rows_to_keep, cols_to_keep]
    }

    output_list <- list(
      propensity_scores = propensity_scores_list,
      summary = ess_summary,
      full_matrix = sorted_matrix
    )

    if (max(study_components$subnetwork) == 2) {
      output_list$subnetwork_matrix <- filtered_matrix
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
        dplyr::group_by(.study) %>%
        dplyr::summarise(
          total_n = dplyr::n(),
          dplyr::across(all_of(covariates), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
          .groups = "drop"
        )
    }

    # Helper Function: Extract AGD Means
    extract_agd_means <- function(agd_df) {
      if (nrow(agd_df) == 0) return(NULL)

      df <- agd_df[, c(".study", ".sample_size"), drop = FALSE]
      retained_covariates <- c()

      for (cov in covariates) {
        distr_obj <- network$integration_code[[cov]]
        args <- distr_obj

        is_binary_dist <- "prob" %in% names(args)

        add_covariate <- TRUE

        if (is_binary_dist) {
          # --- BINARY LOGIC (Bernoulli/Binomial) ---
          prob_col <- as.character(args$prob)

          if (prob_col %in% colnames(agd_df)) {
            p <- agd_df[[prob_col]]
            df[[paste0(cov, "_mean")]] <- p
            # Auto-calculate SD for binary: sqrt(p * (1-p))
            df[[paste0(cov, "_sd")]] <- sqrt(p * (1 - p))
          } else {
            warning(glue::glue("Probability column '{prob_col}' for covariate '{cov}' not found in AgD. Dropped."))
            add_covariate <- FALSE
          }

        } else {

          mean_col <- as.character(args$mean)
          sd_col <- as.character(args$sd)

        # --- Handle Means ---
        if (mean_col %in% colnames(agd_df)) {
          df[[paste0(cov, "_mean")]] <- agd_df[[mean_col]]
        } else if (cov %in% colnames(agd_df)) {
          df[[paste0(cov, "_mean")]] <- agd_df[[cov]]
        } else {
          warning(glue::glue("Mean for covariate '{cov}' not found in AgD. Covariate dropped."))
          add_covariate <- FALSE
        }

        # --- Handle SDs ---
        if (add_covariate) {
          if (sd_col %in% colnames(agd_df)) {
            # Explicit SD column exists -> Use it
            df[[paste0(cov, "_sd")]] <- agd_df[[sd_col]]

          } else {
              # Continuous variable missing SD (or percentage > 1) -> Drop it
              warning(glue::glue("SD column '{sd_col}' for covariate '{cov}' not found in AgD. Dropped (Continuous variable requires explicit SD)."))
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
      stop(paste0(
        "AgD covariate inputs contain missing values:\n",
        paste(" • ", names(lines_by_var), " missing in studies: ", unname(lines_by_var), collapse = "\n"),
        "\nPlease remove these variables from `covariates`"
      ))
    }

    # Weighted Summary of AGD
    agd_summary <- agd_all %>%
      dplyr::group_by(.study) %>%
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

    # 7. Identify Subnetworks
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
      dplyr::select(-.trt) %>%
      dplyr::distinct(.study, subnetwork)

    all_summary <- dplyr::left_join(all_summary, study_components, by = ".study")

    sub1 <- dplyr::filter(all_summary, subnetwork == 1)
    sub2 <- dplyr::filter(all_summary, subnetwork == 2)

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

#' Calculate Latent Bayesian R2 (LOO-Adjusted)
#'
#' Calculates the total Bayesian R2 on the latent scale, automatically handling
#' complex interaction terms (e.g., "age:.trtclass") created by multinma.
#'
#' @param nma A `stan_nma` object.
#'
#' @return A list containing the single LOO-adjusted point estimate.
#' @name cross_validation
#' @export
cross_validation <- function(nma) {

  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("The 'loo' package is required. Please install it with install.packages('loo').")
  }

  if (!inherits(nma, "stan_nma")) {
    stop("Input must be a 'stan_nma' object.")
  }

  # --- Prepare Data ---
  ipd_data <- nma$network$ipd

  # Get the design matrix (X)
  X <- model.matrix(nma$regression, data = ipd_data)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop = FALSE]
  }

  # Get the Posterior Betas
  beta_samples <- as.matrix(nma$stanfit, pars = "beta")
  colnames(beta_samples) <- sub("^beta\\[(.*)\\]$", "\\1", colnames(beta_samples))

  # Initialize Linear Predictor Matrix (Patients x Iterations)
  n_patients <- nrow(X)
  n_iters <- nrow(beta_samples)
  eta_samples <- matrix(0, nrow = n_patients, ncol = n_iters)

  # Determine Interaction Type (Default to 'independent' if NULL)
  int_type <- if (is.null(nma$class_interactions)) "independent" else nma$class_interactions

  # PATH A: Independent Interactions (Fast Matrix Math)
  if (int_type == "independent") {
    common_cols <- intersect(colnames(X), colnames(beta_samples))
    if (length(common_cols) > 0) {
      X_matched <- X[, common_cols, drop = FALSE]
      beta_matched <- beta_samples[, common_cols, drop = FALSE]
      eta_samples <- X_matched %*% t(beta_matched)
    }
  }

  # PATH B: Common/Class Interactions (Robust Loop)
  if (int_type == "common") {
    # Identify Class column
    target_col_name <- if (".trtclass" %in% colnames(ipd_data)) ".trtclass" else ".trt"
    patient_trts <- as.character(ipd_data[[target_col_name]])

    matched_count <- 0

    for (b_name in colnames(beta_samples)) {
      b_vals <- beta_samples[, b_name]
      contribution <- NULL

      if (grepl(":", b_name)) {
        # Interaction Logic
        parts <- strsplit(b_name, ":")[[1]]
        cov_name <- parts[1]
        trt_part <- parts[2]

        if (cov_name %in% colnames(X)) {
          mask <- sapply(patient_trts, function(t) grepl(t, trt_part, fixed = TRUE))
          if (sum(mask, na.rm = TRUE) > 0) {
            contribution <- (X[, cov_name] * mask) %*% t(b_vals)
          }
        }
      } else {
        # Main Effect Logic
        if (b_name %in% colnames(X)) {
          contribution <- X[, b_name] %*% t(b_vals)
        }
      }

      if (!is.null(contribution)) {
        eta_samples <- eta_samples + contribution
        matched_count <- matched_count + 1
      }
    }
  }

  outcome_type <- nma$network$outcome$ipd

  if (outcome_type %in% c("ordered", "binary")) {
    var_res_scalar <- pi^2 / 3
  } else if (outcome_type == "continuous") {
    sigma <- as.matrix(nma$stanfit, pars = "sigma")
    var_res_scalar <- mean(as.vector(sigma^2))
  } else {
    stop(paste("Outcome type", outcome_type, "not supported."))
  }

  # --- 7. LOO-Adjusted R2 ---
  log_lik <- as.matrix(nma$stanfit, pars = "log_lik")
  loo_obj <- suppressWarnings(loo::loo(log_lik, save_psis = TRUE))
  psis_weights <- weights(loo_obj$psis_object, normalize = TRUE, log = FALSE)

  # Calculate Weighted Mean Risk (Eta) for each patient
  loo_eta <- numeric(n_patients)
  for (i in 1:n_patients) {
    loo_eta[i] <- sum(psis_weights[, i] * eta_samples[i, ])
  }

  # Final Formula
  var_fit_loo <- var(loo_eta)
  r2_loo <- var_fit_loo / (var_fit_loo + var_res_scalar)
  r2_percent <- round(r2_loo * 100, 2)

  return(list(r2_percent = r2_percent))
}
