#' The `nma_summary` class
#'
#' The `nma_summary` class contains posterior summary statistics of model
#' parameters or other quantities of interest, and the draws used to obtain
#' these statistics.
#'
#' @rdname nma_summary-class
#' @name nma_summary-class
#' @aliases nma_summary nma_rank_probs
#'
#' @details Objects of class `nma_summary` have the following components:
#'   \describe{
#'   \item{summary}{A data frame containing the computed summary statistics. If
#'   a regression model was fitted with effect modifier interactions with
#'   treatment, these summaries will be study-specific. In this case, the
#'   corresponding study population is indicated in a column named `.study`.}
#'   \item{sims}{A 3D array \[Iteration, Chain, Parameter\] of MCMC
#'   simulations}
#'   \item{studies}{(Optional) A data frame containing study information,
#'   printed along with the corresponding summary statistics if `summary`
#'   contains a `.study` column. Should have a matching `.study` column.}
#'   }
#'
#' The following attributes may also be set:
#'   \describe{
#'   \item{xlab}{Label for x axis in plots, usually either `"Treatment"` or
#'   `"Contrast"`.}
#'   \item{ylab}{Label for y axis in plots, usually used for the scale e.g.
#'   `"log Odds Ratio"`.}
#'   }
#'
#' The subclass `nma_rank_probs` is used by the function
#' [posterior_rank_probs()], and contains posterior rank probabilities. This
#' subclass does not have a `sims` component, as the rank probabilities are
#' themselves posterior summaries of the ranks (i.e. they do not have a
#' posterior distribution). The posterior ranks from which the rank
#' probabilities are calculated may be obtained from [posterior_ranks()].
#'
NULL

#' Methods for `nma_summary` objects
#'
#' The `as.data.frame()`, `as_tibble()`, and `as.tibble()` methods return the
#' posterior summary statistics in a data frame or tibble. The `as.matrix()` and
#' `as.array()` methods return the posterior draws as a matrix or 3D array
#' (Iteration, Chain, Parameter).
#'
#'
#' @param x A `nma_summary` object
#' @param ... Additional arguments passed on to other methods
#' @param digits Integer number of digits to display
#' @param pars Character vector of parameters to display in the printed summary
#' @param include Logical, are parameters named in `pars` included (`TRUE`) or excluded (`FALSE`)
#'
#' @rdname nma_summary-methods
#'
#' @return
#' @export
#'
print.nma_summary <- function(x, ..., digits = 2, pars, include) {
  # Set defaults for pars, include
  if (missing(include)) {
    include <- !missing(pars)
  } else {
    if (!is.logical(include) || length(include) > 1) abort("`include` should be TRUE or FALSE")
  }
  if (missing(pars)) {
    pars <- c("log_lik", "resdev", "fitted",
              "theta_bar_cum", "theta2_bar_cum",
              "mu", "delta", "lp__")
  } else {
    if (!is.character(pars)) abort("`pars` should be a character vector")
  }

  if (!is.numeric(digits) ||
      length(digits) > 1 ||
      trunc(digits) != digits) abort("`digits` must be a single integer")

  x_sum <- tibble::as_tibble(x) %>%
    dplyr::filter(
      if (include) grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
      else !grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
      )

  # Get numeric columns for formatting, handling effective sample sizes and Rhat separately
  num_col <- setdiff(names(x_sum)[purrr::map_lgl(x_sum, is.numeric)],
                     c("n_eff", "Bulk_ESS", "Tail_ESS", "Rhat"))

  x_sum <- dplyr::mutate_at(x_sum, num_col, ~round(., digits))

  if (rlang::has_name(x_sum, "n_eff")) x_sum$n_eff <- round(x_sum$n_eff, 0)
  if (rlang::has_name(x_sum, "Bulk_ESS")) x_sum$Bulk_ESS <- round(x_sum$Bulk_ESS, 0)
  if (rlang::has_name(x_sum, "Tail_ESS")) x_sum$Tail_ESS <- round(x_sum$Tail_ESS, 0)
  if (rlang::has_name(x_sum, "Rhat")) x_sum$Rhat <- round(x_sum$Rhat, max(2, digits))

  x_sum <- tibble::column_to_rownames(x_sum, "parameter")

  # Format summaries nicely by study, if given
  print_study_block <- function(s, info = NULL, ...) {
    this_study <- unique(s$.study)
    sec_header(paste("Study:", this_study))
    cat("\n")
    if (!is.null(info)) {
      cat("Covariate values:\n")
      info %>%
        dplyr::filter(.data$.study == this_study) %>%
        dplyr::select(-.data$.study) %>%
        dplyr::mutate_if(is.numeric, round, digits = digits) %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")
    }
    s %>% dplyr::select(-.data$.study) %>% print(...)
    cat("\n")
  }

  if (rlang::has_name(x_sum, ".study")) {
    by(x_sum, x_sum$.study, print_study_block, info = x$studies, ..., simplify = FALSE)
  } else {
    print(x_sum, ...)
  }
  invisible(x)
}

#' @param geom Character string specifying the `tidybayes` plot geom to use,
#'   default `"pointintervalh"`
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn
#' @rdname nma_summary-methods
#' @export
plot.nma_summary <- function(x, ...,
                             geom = "pointintervalh",
                             ref_line = NA_real_) {
  # Checks
  if (!rlang::is_string(geom))
    abort("`geom` should be a character string specifying the name of a tidybayes geom.")

  tb_geom <- tryCatch(getExportedValue("tidybayes", paste0("stat_", geom)),
    error = function(err) {
      abort(paste("`geom` should be a character string specifying the name of a tidybayes geom:",
                  err, sep = "\n"))
      })


  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  # Is a horizontal geom specified?
  horizontal <- stringr::str_ends(geom, "h")

  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  # Get draws
  draws <- tibble::as_tibble(as.matrix(x))

  if (packageVersion("tidyr") >= "1.0.0") {
    draws <- tidyr::pivot_longer(draws, cols = dplyr::everything(),
                                 names_to = "parameter", values_to = "value")
  } else {
    draws <- tidyr::gather(key = "parameter",
                           value = "value",
                           dplyr::everything())
  }

  if (has_studies <- rlang::has_name(x$summary, ".study")) {
    draws$Study <- forcats::fct_inorder(factor(
      stringr::str_extract(draws$parameter, "(?<=\\[).+(?=\\:)")))
    draws$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(draws$parameter, "(?<=\\: ).+(?=\\])")))
  } else {
    draws$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(draws$parameter, "(?<=\\[).+(?=\\])")))
  }

  if (horizontal) {

    p <- ggplot2::ggplot(draws, ggplot2::aes(y = .data$Treatment, x = .data$value)) +
      ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::scale_y_discrete(p_xlab, limits = rev(levels(draws$Treatment))) +
      ggplot2::xlab(p_ylab)

    if (has_studies) p <- p + ggplot2::facet_grid(Study~.)

  } else {

    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$Treatment, y = .data$value)) +
      ggplot2::geom_hline(yintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::xlab(p_xlab) + ggplot2::ylab(p_ylab)

    if (has_studies) p <- p + ggplot2::facet_grid(.~Study)

  }

  p <- p +
    do.call(tb_geom, args = list(...)) +
    theme_multinma()

  return(p)
}

#' @rdname nma_summary-methods
#' @export
as.data.frame.nma_summary <- function(x, ...) {
  return(as.data.frame(x$summary, ...))
}

#' @rdname nma_summary-methods
#' @importFrom tibble as.tibble
#' @method as.tibble nma_summary
#' @export
as.tibble.nma_summary <- function(x, ...) {
  return(tibble::as.tibble(x$summary, ...))
}

#' @rdname nma_summary-methods
#' @importFrom tibble as_tibble
#' @export
as_tibble.nma_summary <- function(x, ...) {
  return(tibble::as_tibble(x$summary, ...))
}

#' @rdname nma_summary-methods
#' @export
as.array.nma_summary <- function(x, ...) {
  return(x$sims)
}

#' @rdname nma_summary-methods
#' @export
as.matrix.nma_summary <- function(x, ...){
  # Follow approach in rstan:::as.matrix.stanfit
  a <- as.array(x)
  names_a <- dimnames(a)
  dim_a <- dim(a)
  dim(a) <- c(dim_a[1] * dim_a[2], dim_a[3])
  dimnames(a) <- names_a[-2]
  return(a)
}

#' @rdname nma_summary-methods
#' @export
as.array.nma_rank_probs <- function(x, ...) {
  abort(paste("Objects of class `nma_rank_probs` do not contain a 3D MCMC array of simulations, see ?nma_rank_probs.",
              "  - Use as.matrix(), as.data.frame(), or as_tibble() to access the summary rank probabilities",
              "  - Use posterior_ranks() to obtain the posterior ranks themselves", sep = "\n"))
}

#' @rdname nma_summary-methods
#' @export
as.matrix.nma_rank_probs <- function(x, ...){
  df_summary <- tibble::as_tibble(x, ...)
  if (rlang::has_name(df_summary, ".study")) df_summary <- dplyr::select(df_summary, -.data$.study)
  m <- as.matrix(tibble::column_to_rownames(df_summary, "parameter"))
  return(m)
}

#' Compute summary statistics from a 3D MCMC array
#'
#' @param x A 3D MCMC array
#' @param probs Numeric vector of quantiles of interest
#'
#' @return A data frame of computed summaries, one row per parameter
#' @noRd
summary_mcmc_array <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  if (!is.array(x) || length(dim(x)) != 3) abort("Not a 3D MCMC array, [Iterations, Chains, Parameters]")

  pars <- dimnames(x)[[3]]
  p_mean <- apply(x, 3, mean)
  p_sd <- apply(x, 3, sd)
  p_ess_bulk <- apply(x, 3, rstan::ess_bulk)
  p_ess_tail <- apply(x, 3, rstan::ess_tail)
  p_rhat <- apply(x, 3, rstan::Rhat)
  p_se_mean <- p_sd / sqrt(apply(x, 3, rstan:::ess_mean))

  p_quan <- apply(x, 3, quantile, probs = probs)
  if (length(probs) == 1) {
    p_quan <- tibble::tibble(!! paste0(probs*100, "%") := p_quan)
  } else {
    p_quan <- as.data.frame(t(p_quan))
  }

  ss <- tibble::tibble(
    parameter = pars,
    mean = p_mean, se_mean = p_se_mean, sd = p_sd,
    !!! p_quan,
    Bulk_ESS = p_ess_bulk, Tail_ESS = p_ess_tail, Rhat = p_rhat)

  return(ss)
}
