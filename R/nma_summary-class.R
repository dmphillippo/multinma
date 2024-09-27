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
#'   \item{summary}{A data frame containing the computed summary statistics.
#'   Column `.trt` indicates the corresponding treatment, or columns `.trta` and
#'   `.trtb` indicate the corresponding contrast (`.trtb` vs. `.trta`). If a
#'   regression model was fitted with effect modifier interactions with
#'   treatment, these summaries will be study-specific. In this case, the
#'   corresponding study population is indicated in the `.study` column. If a
#'   multinomial model was fitted, the `.category` column indicates the
#'   corresponding category.}
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
#' posterior summary statistics in a data frame or tibble. The `as.matrix()`
#' method returns a matrix of posterior draws. The `as.array()` method returns a
#' 3D array \[Iteration, Chain, Parameter\] of posterior draws (as class
#' [mcmc_array]).
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
#' @seealso [plot.nma_summary()]
#'
#' @return A `data.frame` for `as.data.frame()`, a `tbl_df` for `as.tibble()`
#'   and `as_tibble()`, a `matrix` for `as.matrix()`, and an `mcmc_array` for
#'   `as.array()`.
#'
#'   The `print()` method returns `x` invisibly.
#' @export
#'
print.nma_summary <- function(x, ..., digits = 2, pars, include = TRUE) {
  # Checks
  if (!rlang::is_bool(include)) abort("`include` should be TRUE or FALSE")

  if (!missing(pars) && !is.character(pars)) abort("`pars` should be a character vector")

  if (!rlang::is_scalar_integerish(digits)) abort("`digits` must be a single integer")

  x_sum <- tibble::as_tibble(x)

  if (!missing(pars)) {
    x_sum <-
      dplyr::filter(x_sum,
        if (include) grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
        else !grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
        )
  }

  # Get numeric columns for formatting, handling effective sample sizes and Rhat separately
  num_col <- setdiff(names(x_sum)[purrr::map_lgl(x_sum, is.numeric)],
                     c("n_eff", "Bulk_ESS", "Tail_ESS", "Rhat"))

  x_sum <- dplyr::mutate_at(x_sum, num_col, ~round(., digits))

  if (rlang::has_name(x_sum, "n_eff")) x_sum$n_eff <- round(x_sum$n_eff, 0)
  if (rlang::has_name(x_sum, "Bulk_ESS")) x_sum$Bulk_ESS <- round(x_sum$Bulk_ESS, 0)
  if (rlang::has_name(x_sum, "Tail_ESS")) x_sum$Tail_ESS <- round(x_sum$Tail_ESS, 0)
  if (rlang::has_name(x_sum, "Rhat")) x_sum$Rhat <- round(x_sum$Rhat, max(2, digits))

  x_sum <- tibble::column_to_rownames(x_sum, "parameter")

  # Drop dot columns (.trta, .trtb, .trt, .category) from the output, except .study and .time
  x_sum <- dplyr::select(x_sum, -dplyr::matches("^\\.(?!study|time)", perl = TRUE))

  # Format summaries nicely by study, if given
  print_study_block <- function(s, info = NULL, ...) {
    this_study <- unique(s$.study)
    sec_header(paste("Study:", this_study))
    cat("\n")
    if (!is.null(info)) {
      cat("Covariate values:\n")
      info %>%
        dplyr::filter(.data$.study == this_study) %>%
        dplyr::select(-".study") %>%
        dplyr::mutate_if(is.numeric, round, digits = digits) %>%
        as.data.frame() %>%
        print(..., row.names = FALSE)
      cat("\n")
    }
    s %>% dplyr::select(-".study") %>% print(...)
    cat("\n")
  }

  if (rlang::has_name(x_sum, ".study")) {
    by(x_sum, x_sum$.study, print_study_block, info = x$studies, ..., simplify = FALSE)
  } else {
    print(x_sum, ...)
  }
  invisible(x)
}

#' Plots of summary results
#'
#' The `plot` method for `nma_summary` objects is used to produce plots of
#' parameter estimates (when called on a `stan_nma` object or its summary),
#' relative effects (when called on the output of [relative_effects()]),
#' absolute predictions (when called on the output of [predict.stan_nma()]),
#' posterior ranks and rank probabilities (when called on the output of
#' [posterior_ranks()] or [posterior_rank_probs()]).
#'
#' @param x A `nma_summary` object
#' @param ... Additional arguments passed on to the underlying `ggdist` plot
#'   stat, see Details
#' @param stat Character string specifying the `ggdist` plot stat to use,
#'   default `"pointinterval"`, except when plotting estimated
#'   survival/hazard/cumulative hazard curves from survival models where the
#'   default is `"lineribbon"`
#' @param orientation Whether the `ggdist` geom is drawn horizontally
#'   (`"horizontal"`) or vertically (`"vertical"`), default `"horizontal"`
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn
#'
#' @details Plotting is handled by \link[ggplot2:ggplot2-package]{ggplot2} and
#'   the stats and geoms provided in the \link[ggdist:ggdist-package]{ggdist}
#'   package. As a result, the output is very flexible. Any plotting stats
#'   provided by `ggdist` may be used, via the argument `stat`.
#'
#'   The default uses
#'   \code{\link[ggdist:stat_pointinterval]{ggdist::stat_pointinterval()}}, to
#'   produce medians and 95% Credible Intervals with 66% inner bands. Additional
#'   arguments in `...` are passed to the `ggdist` stat, to customise the
#'   output. For example, to produce means and Credible Intervals, specify
#'   `point_interval = "mean_qi"`. To produce an 80% Credible Interval with no
#'   inner band, specify `.width = c(0, 0.8)`.
#'
#'   Alternative stats can be specified to produce different summaries. For
#'   example, specify `stat = "[half]eye"` to produce (half) eye plots, or `stat
#'   = "histinterval"` to produce histograms with intervals.
#'
#'   A full list of options and examples is found in the `ggdist` vignette
#'   `vignette("slabinterval", package = "ggdist")`.
#'
#'   For survival/hazard/cumulative hazard curves estimated from survival
#'   models, the default uses
#'   \code{\link[ggdist:stat_lineribbon]{ggdist::stat_lineribbon()}} which
#'   produces curves of posterior medians with 50%, 80%, and 95% Credible
#'   Interval bands. Again, additional arguments in `...` are passed to the
#'   `ggdist` stat. For example, to produce posterior means and 95% Credible
#'   bands, specify `point_interval = "mean_qi"` and `.width = 0.95`.
#'
#'   A `ggplot` object is returned which can be further modified through the
#'   usual \link[ggplot2:ggplot2-package]{ggplot2} functions to add further
#'   aesthetics, geoms, themes, etc.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @examples
#' ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Produce relative effects
#' smk_releff_RE <- relative_effects(smk_fit_RE)
#' plot(smk_releff_RE, ref_line = 0)
#'
#' # Customise plot options
#' plot(smk_releff_RE, ref_line = 0, stat = "halfeye")
#'
#' # Further customisation is possible with ggplot commands
#' plot(smk_releff_RE, ref_line = 0, stat = "halfeye", slab_alpha = 0.6) +
#'   ggplot2::aes(slab_fill = ggplot2::after_stat(ifelse(x < 0, "darkred", "grey60")))
#'
#' # Produce posterior ranks
#' smk_rank_RE <- posterior_ranks(smk_fit_RE, lower_better = FALSE)
#' plot(smk_rank_RE)
#'
#' # Produce rank probabilities
#' smk_rankprob_RE <- posterior_rank_probs(smk_fit_RE, lower_better = FALSE)
#' plot(smk_rankprob_RE)
#'
#' # Produce cumulative rank probabilities
#' smk_cumrankprob_RE <- posterior_rank_probs(smk_fit_RE, lower_better = FALSE,
#'                                            cumulative = TRUE)
#' plot(smk_cumrankprob_RE)
#'
#' # Further customisation is possible with ggplot commands
#' plot(smk_cumrankprob_RE) +
#'   ggplot2::facet_null() +
#'   ggplot2::aes(colour = Treatment)
#' }
plot.nma_summary <- function(x, ...,
                             stat = "pointinterval",
                             orientation = c("horizontal", "vertical", "y", "x"),
                             ref_line = NA_real_) {
  # Checks
  if (!rlang::is_string(stat))
    abort("`stat` should be a character string specifying the name of a ggdist stat")

  stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

  tb_geom <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
    error = function(err) {
      abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                  err, sep = "\n"))
      })


  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  orientation <- rlang::arg_match(orientation)
  if (orientation == "x") orientation <- "vertical"
  else if (orientation == "y") orientation <- "horizontal"

  # Is a horizontal geom specified?
  horizontal <- orientation == "horizontal"

  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  # Tweak output if rank summary
  if (is_ranks <- inherits(x, "nma_ranks")) {
    ntrt <- nrow(as.data.frame(x))
  }

  # Get draws
  draws <- tibble::as_tibble(as.matrix(x))

  draws <- tidyr::pivot_longer(draws, cols = dplyr::everything(),
                               names_to = "parameter", values_to = "value")

  if (has_studies <- rlang::has_name(as.data.frame(x), ".study")) {
    draws$Study <- forcats::fct_inorder(factor(
      stringr::str_extract(draws$parameter, "(?<=\\[).+(?=\\:)")))

    if (inherits(x, "ordered_nma_summary")) {
      draws$Treatment <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=\\: ).+(?=, .+?\\])")))
      draws$Category <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=, ).+(?=\\])")))
    } else {
      draws$Treatment <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=\\: ).+(?=\\])")))
    }
  } else {
    if (inherits(x, "ordered_nma_summary")) {
      draws$Treatment <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=\\[).+(?=, .+?\\])")))
      draws$Category <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=, ).+(?=\\])")))
    } else {
      draws$Treatment <- forcats::fct_inorder(factor(
        stringr::str_extract(draws$parameter, "(?<=\\[).+(?=\\])")))
    }
  }

  if (horizontal) {

    draws$Treatment <- forcats::fct_rev(draws$Treatment)

    p <- ggplot2::ggplot(draws, ggplot2::aes(y = .data$Treatment, x = .data$value)) +
      ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::ylab(p_xlab)

    if (has_studies) {
      if (inherits(x, "ordered_nma_summary")) {
        p <- p + ggplot2::facet_grid(Study~Category)
      } else {
        p <- p + ggplot2::facet_grid(Study~.)
      }
    } else if (inherits(x, "ordered_nma_summary")) {
      p <- p + ggplot2::facet_grid(.~Category)
    }

    if (is_ranks) {
      p <- p + ggplot2::scale_x_continuous(p_ylab, breaks = 1:ntrt, minor_breaks = NULL)
    } else {
      p <- p + ggplot2::xlab(p_ylab)
    }

  } else {

    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$Treatment, y = .data$value)) +
      ggplot2::geom_hline(yintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::xlab(p_xlab)

    if (has_studies) {
      if (inherits(x, "ordered_nma_summary")) {
        p <- p + ggplot2::facet_grid(Category~Study)
      } else {
        p <- p + ggplot2::facet_grid(.~Study)
      }
    } else if (inherits(x, "ordered_nma_summary")) {
      p <- p + ggplot2::facet_grid(Category~.)
    }


    if (is_ranks) {
      p <- p + ggplot2::scale_y_continuous(p_ylab, breaks = 1:ntrt, minor_breaks = NULL)
    } else {
      p <- p + ggplot2::ylab(p_ylab)
    }

  }

  p <- p +
    do.call(tb_geom, args = list(orientation = orientation, ...)) +
    theme_multinma()

  return(p)
}

#' @rdname plot.nma_summary
#' @export
plot.nma_parameter_summary <- function(x, ...,
                                       stat = "pointinterval",
                                       orientation = c("horizontal", "vertical", "y", "x"),
                                       ref_line = NA_real_) {
  # Checks
  if (!rlang::is_string(stat))
    abort("`stat` should be a character string specifying the name of a ggdist stat.")

  stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

  tb_geom <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
                      error = function(err) {
                        abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                                    err, sep = "\n"))
                      })


  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  orientation <- rlang::arg_match(orientation)
  if (orientation == "x") orientation <- "vertical"
  else if (orientation == "y") orientation <- "horizontal"

  # Is a horizontal geom specified?
  horizontal <- orientation == "horizontal"

  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  # Get draws
  draws <- tibble::as_tibble(as.matrix(x))

  draws <- tidyr::pivot_longer(draws, cols = dplyr::everything(),
                               names_to = "parameter", values_to = "value")

  draws$par_base <- forcats::fct_inorder(factor(
    stringr::str_remove(draws$parameter, "\\[.*\\]")))
  draws$parameter <- forcats::fct_inorder(factor(draws$parameter))

  if (horizontal) {

    draws$parameter <- forcats::fct_rev(draws$parameter)

    p <- ggplot2::ggplot(draws, ggplot2::aes(y = .data$parameter, x = .data$value)) +
      ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::ylab(p_xlab) + ggplot2::xlab(p_ylab) +
      ggplot2::facet_grid(par_base~., scales = "free", space = "free")

  } else {

    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$parameter, y = .data$value)) +
      ggplot2::geom_hline(yintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::xlab(p_xlab) + ggplot2::ylab(p_ylab) +
      ggplot2::facet_grid(.~par_base, scales = "free", space = "free")

  }

  p <- p +
    do.call(tb_geom, args = list(orientation = orientation, ...)) +
    theme_multinma()

  return(p)
}

#' @rdname plot.nma_summary
#' @export
plot.nma_rank_probs <- function(x, ...) {
  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  dat <- as.data.frame(x)

  ntrt <- nrow(dat)

  if (has_studies <- rlang::has_name(dat, ".study")) {
    dat$Study <- forcats::fct_inorder(factor(dat$.study))
    dat$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(dat$parameter, "(?<=\\: ).+(?=\\])")))
  } else {
    dat$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(dat$parameter, "(?<=\\[).+(?=\\])")))
  }

  dat <- tidyr::pivot_longer(dat, cols = dplyr::starts_with("p_rank"),
                               names_to = "rank", values_to = "probability",
                               names_pattern = "^p_rank\\[([0-9]+)\\]$",
                               names_transform = list(rank = as.integer))


  p <- ggplot2::ggplot(dat,
                       ggplot2::aes(x = .data$rank, y = .data$probability)) +
    ggplot2::geom_line(...) +
    ggplot2::ylab(p_ylab) +
    ggplot2::scale_x_continuous(p_xlab, breaks = 1:ntrt, minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_multinma()

  if (has_studies) {
    p <- p + ggplot2::facet_grid(Study~Treatment)
  } else {
    p <- p + ggplot2::facet_wrap(~Treatment)
  }

  return(p)
}

#' @rdname plot.nma_summary
#' @export
plot.surv_nma_summary <- function(x, ..., stat = "lineribbon") {
  # Checks
  if (!rlang::is_string(stat))
    abort("`stat` should be a character string specifying the name of a ggdist stat")

  stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

  tb_geom <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
                      error = function(err) {
                        abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                                    err, sep = "\n"))
                      })

  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  # Get draws
  draws <- as.data.frame(x) %>%
    dplyr::select(dplyr::starts_with(".")) %>%
    dplyr::mutate(value = purrr::array_branch(as.matrix(x), 2)) %>%
    tidyr::unnest(cols = "value")

  if (has_studies <- rlang::has_name(draws, ".study")) {
    draws$Study <- draws$.study
  }

  draws$Time <- draws$.time

  if (rlang::has_name(draws, ".trt") || dplyr::n_distinct(draws$.trta) == 1) {
    draws$Treatment <- if (rlang::has_name(draws, ".trt")) draws$.trt else draws$.trtb
    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$Time, y = .data$value,
                                             colour = .data$Treatment,
                                             fill = .data$Treatment))
  } else {
    draws$Contrast <- paste0(draws$.trtb, " vs. ", draws$.trta)
    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$Time, y = .data$value,
                                             colour = .data$Contrast,
                                             fill = .data$Contrast))
  }

  p <- p +
    ggplot2::ylab(p_ylab) +
    ggplot2::xlab(p_xlab)

  if (has_studies) {
      p <- p + ggplot2::facet_grid(.~Study)
  }

  p <- p +
    # Draw lines and intervals separately for better clarity
    do.call(tb_geom, args = rlang::dots_list(linewidth = 0, ..., alpha = 0.15, .homonyms = "first")) +
    do.call(tb_geom, args = rlang::dots_list(.width = 0, alpha = 1, ..., linewidth = 0.5, .homonyms = "first")) +
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
  out <- x$sims
  class(out) <-  c("mcmc_array", "array")
  return(out)
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
  class(a) <- "matrix"
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
  if (rlang::has_name(df_summary, ".study")) df_summary <- dplyr::select(df_summary, -".study")
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
  check_probs(probs)

  pars <- dimnames(x)[[3]]
  p_mean <- apply(x, 3, mean)
  p_sd <- apply(x, 3, sd)
  p_ess_bulk <- apply(x, 3, rstan::ess_bulk)
  p_ess_tail <- apply(x, 3, rstan::ess_tail)
  p_rhat <- apply(x, 3, rstan::Rhat)
  # p_se_mean <- p_sd / sqrt(apply(x, 3, rstan:::ess_mean))

  qt <- function(x, probs, ...) {
    if (all(is.na(x))) setNames(rlang::rep_along(probs, NA_real_), paste0(probs*100, "%"))
    else quantile(x, probs = probs, ...)
  }

  p_quan <- apply(x, 3, qt, probs = probs)
  if (length(probs) == 1) {
    p_quan <- tibble::tibble(!! paste0(probs*100, "%") := p_quan)
  } else {
    p_quan <- as.data.frame(t(p_quan))
  }

  ss <- tibble::tibble(
    parameter = pars,
    mean = p_mean,
    # se_mean = p_se_mean,
    sd = p_sd,
    !!! p_quan,
    Bulk_ESS = p_ess_bulk, Tail_ESS = p_ess_tail, Rhat = p_rhat)

  return(ss)
}

#' Validate probs argument
#' @noRd
check_probs <- function(probs) {
  if (!rlang::is_double(probs, finite = TRUE) || any(probs < 0) || any(probs > 1))
    rlang::abort("`probs` must be a numeric vector of probabilities between 0 and 1.")
}
