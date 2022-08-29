#' Treatment rankings and rank probabilities
#'
#' Produce posterior treatment rankings and rank probabilities from a fitted NMA
#' model. When a meta-regression is fitted with effect modifier interactions
#' with treatment, these will differ by study population.
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
#' @param sucra Logical, calculate the surface under the cumulative ranking
#'   curve (SUCRA) for each treatment? Default `FALSE`.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @seealso [plot.nma_summary()] for plotting the ranks and rank probabilities.
#'
#' @details The function `posterior_ranks()` produces posterior rankings, which
#'   have a distribution (e.g. mean/median rank and 95% Credible Interval). The
#'   function `posterior_rank_probs()` produces rank probabilities, which give
#'   the posterior probabilities of being ranked first, second, etc. out of all
#'   treatments.
#'
#'   The argument `lower_better` specifies whether lower treatment
#'   effects or higher treatment effects are preferred. For example, with a
#'   negative binary outcome lower (more negative) log odds ratios are
#'   preferred, so `lower_better = TRUE`. Conversely, for example, if treatments
#'   aim to increase the rate of a positive outcome then `lower_better = FALSE`.
#'
#' @examples
#' ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Produce posterior ranks
#' smk_rank_RE <- posterior_ranks(smk_fit_RE, lower_better = FALSE)
#' smk_rank_RE
#' plot(smk_rank_RE)
#'
#' # Produce rank probabilities
#' smk_rankprob_RE <- posterior_rank_probs(smk_fit_RE, lower_better = FALSE)
#' smk_rankprob_RE
#' plot(smk_rankprob_RE)
#'
#' # Produce cumulative rank probabilities
#' smk_cumrankprob_RE <- posterior_rank_probs(smk_fit_RE, lower_better = FALSE,
#'                                            cumulative = TRUE)
#' smk_cumrankprob_RE
#' plot(smk_cumrankprob_RE)
#'
#' # Further customisation is possible with ggplot commands
#' plot(smk_cumrankprob_RE) +
#'   ggplot2::facet_null() +
#'   ggplot2::aes(colour = Treatment)
#' }
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_mlnmr_example
#' @examples \donttest{
#' # Produce population-adjusted rankings for all study populations in
#' # the network
#'
#' # Ranks
#' pso_rank <- posterior_ranks(pso_fit)
#' pso_rank
#' plot(pso_rank)
#'
#' # Rank probabilities
#' pso_rankprobs <- posterior_rank_probs(pso_fit)
#' pso_rankprobs
#' plot(pso_rankprobs)
#'
#' # Cumulative rank probabilities
#' pso_cumrankprobs <- posterior_rank_probs(pso_fit, cumulative = TRUE)
#' pso_cumrankprobs
#' plot(pso_cumrankprobs)
#'
#' # Produce population-adjusted rankings for a different target
#' # population
#' new_agd_means <- data.frame(
#'   bsa = 0.6,
#'   prevsys = 0.1,
#'   psa = 0.2,
#'   weight = 10,
#'   durnpso = 3)
#'
#' # Ranks
#' posterior_ranks(pso_fit, newdata = new_agd_means)
#'
#' # Rank probabilities
#' posterior_rank_probs(pso_fit, newdata = new_agd_means)
#'
#' # Cumulative rank probabilities
#' posterior_rank_probs(pso_fit, newdata = new_agd_means,
#'                      cumulative = TRUE)
#' }
posterior_ranks <- function(x, newdata = NULL, study = NULL,
                            lower_better = TRUE,
                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            sucra = FALSE,
                            summary = TRUE) {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")

  check_probs(probs)

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce ranks under inconsistency '{x$consistency}' model."))

  # Get reference treatment, number of treatments
  trt_ref <- levels(x$network$treatments)[1]
  ntrt <- nlevels(x$network$treatments)

  # All other checks handled by relative_effects()
  rel_eff <- relative_effects(x = x, newdata = newdata, study = {{ study }},
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
    rk <- aperm(apply(d, 1:2, rank, ties.method = "min"),
                c("iterations", "chains", "parameters"))

    # If higher treatment effects are better
    if (!lower_better) rk <- ntrt + 1 - rk

    # Rename parameters
    dimnames(rk)[[3]] <- paste0("rank[", levels(x$network$treatments), "]")

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs) %>%
        tibble::add_column(.trt = x$network$treatments, .before = 1)

      if (sucra) {
        # Calculate SUCRA using scaled mean rank relation of Rucker and Schwarzer (2015)
        sucras <- unname((ntrt - rk_summary$mean) / (ntrt - 1))
        rk_summary <- tibble::add_column(rk_summary, sucra = sucras, .after = "sd")
      }

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
        aperm(apply(temp_d, 1:2, rank, ties.method = "min"),
              c("iterations", "chains", "parameters"))
    }

    # If higher treatment effects are better
    if (!lower_better) rk <- ntrt + 1 - rk

    # Rename parameters
    dimnames(rk)[[3]] <- gsub("^d\\[", "rank\\[", rk_names)

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs) %>%
        # Add in study info
        tibble::add_column(.study = rep(studies$.study, each = ntrt),
                           .trt = rep(x$network$treatments, times = nstudy),
                           .before = 1)

      if (sucra) {
        # Calculate SUCRA using scaled mean rank relation of Rucker and Schwarzer (2015)
        sucras <- unname((ntrt - rk_summary$mean) / (ntrt - 1))
        rk_summary <- tibble::add_column(rk_summary, sucra = sucras, .after = "sd")
      }

      out <- list(summary = rk_summary, sims = rk, studies = studies)
    } else {
      out <- list(sims = rk, studies = studies)
    }
  }

  if (summary) {
    class(out) <- c("nma_ranks", "nma_summary")
    attr(out, "xlab") <- "Treatment"
    attr(out, "ylab") <- "Posterior Rank"
  }
  return(out)
}

#' @param cumulative Logical, return cumulative rank probabilities? Default is
#'   `FALSE`, return posterior probabilities of each treatment having a given
#'   rank. If `TRUE`, cumulative posterior rank probabilities are returned for
#'   each treatment having a given rank or better.
#' @export
#' @rdname posterior_ranks
posterior_rank_probs <- function(x, newdata = NULL, study = NULL, lower_better = TRUE,
                                 cumulative = FALSE, sucra = FALSE) {
  # Checks
  if (!rlang::is_bool(cumulative))
    abort("`cumulative` should be TRUE or FALSE.")

  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")

  # All other checks handled by posterior_ranks()
  rk <- posterior_ranks(x = x, newdata = newdata, study = {{ study }},
                        lower_better = lower_better, summary = FALSE)

  ntrt <- nlevels(x$network$treatments)
  studies <- rk$studies

  if (is.null(studies)) { # No study-specific treatment effects

    p_rank <- apply(apply(rk$sims, 1:3, `==`, 1:ntrt), c(4, 1), mean)
    if (cumulative) p_rank <- t(apply(p_rank, 1, cumsum))

    if (sucra) {
      if (cumulative) sucras <- rowMeans(p_rank[, -ntrt, drop = FALSE])
      else sucras <- rowMeans(t(apply(p_rank, 1, cumsum))[, -ntrt, drop = FALSE])
    }

    rownames(p_rank) <- stringr::str_replace(rownames(p_rank), "^rank\\[", "d\\[")
    colnames(p_rank) <- paste0("p_rank[", 1:ntrt, "]")

    p_rank <- tibble::as_tibble(p_rank, rownames = "parameter") %>%
      tibble::add_column(.trt = x$network$treatments, .before = 1)

    if (sucra) p_rank$sucra <- unname(sucras)

    out <- list(summary = p_rank)

  } else { # Study-specific treatment effects

    nstudy <- nrow(studies)

    p_rank <- matrix(NA_real_, nrow = ntrt * nstudy, ncol = ntrt)
    rk_sims <- rk$sims
    for (i in seq_len(nstudy)) {
      p_rank[(i - 1)*ntrt + 1:ntrt, ] <-
        apply(apply(rk_sims[ , , (i - 1)*ntrt + 1:ntrt], 1:3, `==`, 1:ntrt), c(4, 1), mean)
    }

    if (cumulative) p_rank <- t(apply(p_rank, 1, cumsum))

    if (sucra) {
      if (cumulative) sucras <- rowMeans(p_rank[, -ntrt, drop = FALSE])
      else sucras <- rowMeans(t(apply(p_rank, 1, cumsum))[, -ntrt, drop = FALSE])
    }

    rownames(p_rank) <- stringr::str_replace(dimnames(rk_sims)[[3]], "^rank\\[", "d\\[")
    colnames(p_rank) <- paste0("p_rank[", 1:ntrt, "]")

    p_rank <- tibble::as_tibble(p_rank, rownames = "parameter") %>%
      # Add in study info
      tibble::add_column(.study = rep(studies$.study, each = ntrt),
                         .trt = rep(x$network$treatments, times = nstudy),
                         .before = 1)

    if (sucra) p_rank$sucra <- unname(sucras)

    out <- list(summary = p_rank, studies = studies)

  }

  class(out) <- c("nma_rank_probs", "nma_summary")
  attr(out, "xlab") <- "Rank"
  attr(out, "ylab") <- if (cumulative) "Cumulative Rank Probability" else "Rank Probability"
  return(out)
}
