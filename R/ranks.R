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
#'
#' @return A [nma_summary] object
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
                           probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce ranks under inconsistency '{x$consistency}' model."))

  # Get reference treatment, number of treatments
  trt_ref <- levels(x$network$treatments)[1]
  ntrt <- nlevels(x$network$treatments)

  # All other checks handled by relative_effects()
  rel_eff <- relative_effects(x = x, newdata = newdata, study = enquo(study),
                              all_contrasts = FALSE, summary = FALSE)


  if (is.null(rel_eff$studies)) { # No study-specific treatment effects
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
    rk_summary <- summary_mcmc_array(rk, probs = probs)
    out <- list(summary = rk_summary, sims = rk)

  } else { # Study-specific treatment effects

  }

  class(out) <- "nma_summary"
  return(out)
}

#' @export
#' @rdname posterior_rank
#' @examples
posterior_rank_prob <- function(x, newdata = NULL, study = NULL, lower_better = TRUE) {

}
