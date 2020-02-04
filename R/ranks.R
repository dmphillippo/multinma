#' Treatment rankings
#'
#' @param x A `stan_nma` object created by [nma()]
#' @param newdata Only used if a regression model is fitted. A data frame of
#'   study details, one row per study, giving the covariate values at which to
#'   produce relative effects. Column names must match variables in the
#'   regression model. If `NULL`, relative effects are produced for all studies
#'   in the network.
#' @param study Column of `newdata` which specifies study names, otherwise
#'   studies will be labelled by row number.
#' @param lower_better Logical, are lower treatment effects better (`TRUE`;
#'   default) or higher better (`FALSE`)? See details.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @details The argument `lower_better` specifies whether lower treatment
#'   effects or higher treatment effects are preferred. For example, with a
#'   negative binary outcome lower (more negative) log odds ratios are
#'   preferred, so `lower_better = TRUE`. Conversely, for example, if treatments
#'   aim to increase the rate of a positive outcome then `lower_better =
#'   FALSE`.
#'
#' @examples
posterior_rank <- function(x, newdata = NULL, study = NULL,
                           lower_better = TRUE,
                           probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                           summary = TRUE) {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce ranks under inconsistency '{x$consistency}' model."))

  # Get reference treatment, number of treatments
  trt_ref <- levels(x$network$treatments)[1]
  ntrt <- nlevels(x$network$treatments)

  # All other checks handled by relative_effects()
  rel_eff <- relative_effects(x = x, newdata = newdata, study = enquo(study),
                              all_contrasts = FALSE, summary = FALSE)

  studies <- rel_eff$studies

  if (is.null(studies)) { # No study-specific treatment effects
    # Add zeros for d[1]
    dim_d <- dim(rel_eff$sim)
    dim_d[3] <- dim_d[3] + 1
    dimnames_d <- dimnames(rel_eff$sim)
    dimnames_d[[3]] <- c(paste0("d[", trt_ref, "]"), dimnames_d[[3]])
    d <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
    d[ , , 1] <- 0
    d[ , , 2:ntrt] <- rel_eff$sim

    # Get ranks at each iteration
    rk <- aperm(apply(d, 1:2, rank, ties.method = "average"),
                c("iterations", "chains", "parameters"))

    # Rename parameters
    dimnames(rk)[[3]] <- paste0("rank[", levels(x$network$treatments), "]")

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs)
      out <- list(summary = rk_summary, sims = rk)
    } else {
      out <- list(sims = rk)
    }

  } else { # Study-specific treatment effects
    nstudy <- nrow(studies)

    d <- rel_eff$sim
    d_names <- dimnames(d)[[3]]

    # Calculate ranks within each study population
    dim_rk <- dim(d)
    dim_rk[3] <- dim_rk[3] + nstudy
    rk_names <- vector("character", dim_rk[[3]])
    dimnames_rk <- dimnames(d)
    dimnames_rk[[3]] <- rk_names
    rk <- array(NA_real_, dim = dim_rk, dimnames = dimnames_rk)

    dim_temp_d <- dim(d)
    dim_temp_d[[3]] <- ntrt
    dimnames_temp_d <- dimnames(d)
    dimnames_temp_d[[3]] <- paste0("d[", levels(x$network$treatments), "]")
    temp_d <- array(NA_real_, dim = dim_temp_d, dimnames = dimnames_temp_d)
    temp_d[ , , 1] <- 0

    for (i in seq_len(nstudy)) {
      rk_names[(i - 1)*ntrt + 1:ntrt] <-
        c(paste0("d[", studies$.study[i], ": ", trt_ref, "]"),
          d_names[(i - 1)*(ntrt - 1) + 1:(ntrt - 1)])

      temp_d[ , , 2:ntrt] <- d[ , , (i - 1)*(ntrt - 1) + 1:(ntrt - 1)]

      rk[ , , (i - 1)*ntrt + 1:ntrt] <-
        aperm(apply(temp_d, 1:2, rank, ties.method = "average"),
              c("iterations", "chains", "parameters"))
    }

    # Rename parameters
    dimnames(rk)[[3]] <- gsub("^d\\[", "rank\\[", rk_names)

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs) %>%
        # Add in study info
        tibble::add_column(.study = rep(studies$.study, each = ntrt), .before = 1)

      out <- list(summary = rk_summary, sims = rk, studies = studies)
    } else {
      out <- list(sims = rk, studies = studies)
    }
  }

  if (summary) class(out) <- "nma_summary"
  return(out)
}

#' @export
#' @rdname posterior_rank
#' @examples
posterior_rank_prob <- function(x, newdata = NULL, study = NULL, lower_better = TRUE) {

}
