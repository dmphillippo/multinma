#' The stan_nma class
#'
#' The `stan_nma` and `stan_mlnmr` classes contains the results from running a
#' model with the function [nma()].
#'
#' @rdname stan_nma-class
#' @name stan_nma-class
#' @aliases stan_nma stan_mlnmr
#'
#' @details Objects of class `stan_nma` and `stan_mlnmr` have the following
#'   components:
#'   \describe{
#'   \item{`network`}{The network data from which the model was run (class
#'   [nma_data] for `stan_nma`, or class [mlnmr_data] for `stan_mlnmr`)}
#'   \item{`stanfit`}{The `stanfit` object returned by calling
#'   `\link[rstan:stanmodel-method-sampling]{sampling()}` for the model}
#'   \item{`trt_effects`}{Whether fixed or random effects were used (character
#'   string)}
#'   \item{`consistency`}{The consistency/inconsistency model used (character
#'   string)}
#'   \item{`regression`}{The regression model used (formula)}
#'   \item{`class_interactions`}{If treatment classes and a regression model are
#'   specified, the model used for interactions within each class (common,
#'   exchangeable, or independent)}
#'   \item{`xbar`}{A named vector of values used for centering}
#'   \item{`likelihood`}{The likelihood used (character string)}
#'   \item{`link`}{The link function used (character string)}
#'   \item{`priors`}{A list containing the priors used (as [nma_prior] objects)}
#'   }
#'
#' The `stan_mlnmr` sub-class inherits from `stan_nma`, and differs only in the
#' class of the `network` object.
#'
NULL

#' Print `stan_nma` objects
#'
#' @param x A [stan_nma] object
#' @param ... Further arguments passed to [print.stanfit()]
#'
#' @export
print.stan_nma <- function(x, ...) {
  if (inherits(x$network, "mlnmr_data")) type <- "ML-NMR"
  else type <- "NMA"
  cglue("A {x$trt_effects} effects {type} with a {x$likelihood} likelihood ({x$link} link).")
  if (x$consistency != "consistency") cglue("An inconsistency model ('{x$consistency}') was fitted.")
  if (!is.null(x$regression)) {
    cglue("Regression model: {rlang::as_label(x$regression)}.")
    if (!is.null(x$xbar)) {
      cglue("Centred covariates at the following overall mean values:")
      print(x$xbar)
    }
  }

  sf <- as.stanfit(x)
  dots <- list(...)
  include <- "pars" %in% names(dots)
  dots <- rlang::dots_list(x = sf,
                           pars = c("log_lik", "resdev", "fitted",
                                    "theta_bar_cum", "theta2_bar_cum",
                                    "mu", "delta"),
                           include = include,
                           use_cache = FALSE,
                           !!! dots,
                           .homonyms = "last")
  do.call(print, dots)
  invisible(x)
}

#' Posterior summaries from `stan_nma` objects
#'
#' Posterior summaries of model parameters in `stan_nma` objects may be produced
#' using the `summary()` method and plotted with the `plot()` method. NOTE: To
#' produce relative effects, absolute predictions, or posterior ranks, see
#' [relative_effects()], [predict.stan_nma()], [posterior_ranks()],
#' [posterior_rank_probs()].
#'
#' @param x A `stan_nma` object
#' @param ... Additional arguments passed on to other methods
#' @param pars,include See `\link[rstan:stanfit-method-extract]{rstan::extract()}`
#' @param probs Numeric vector of specifying quantiles of interest, default
#'   `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#'
#' @details The `plot()` method is a shortcut for `plot(summary(stan_nma))`. For
#'   details of plotting options, see [plot.nma_summary()].
#'
#' @return A [nma_summary] object
#' @export
#'
#' @seealso [plot.nma_summary()], [relative_effects()], [predict.stan_nma()],
#'   [posterior_ranks()], [posterior_rank_probs()]
#'
#' @examples
summary.stan_nma <- function(x, ...,
                             pars, include,
                             probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
                             ) {

  # Set defaults for pars, include
  if (missing(include)) {
    include <- !missing(pars)
  } else {
    if (!rlang::is_bool(include)) abort("`include` should be TRUE or FALSE")
  }
  if (missing(pars)) {
    pars <- c("log_lik", "resdev", "fitted",
              "theta_bar_cum", "theta2_bar_cum",
              "lp__")
  } else {
    if (!is.character(pars)) abort("`pars` should be a character vector")
  }

  sims <- as.array(x, pars = pars, include = include)
  sums <- summary_mcmc_array(sims, probs = probs)
  ss <- list(summary = sums, sims = sims)
  class(ss) <- c("nma_parameter_summary", "nma_summary")
  attr(ss, "xlab") <- "Parameter"
  attr(ss, "ylab") <- "Value"
  return(ss)
}

#' @param stat Character string specifying the `tidybayes` plot stat to use,
#'   default `"pointintervalh"`
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn
#' @rdname summary.stan_nma
#' @export
plot.stan_nma <- function(x, ...,
                          pars, include,
                          stat = "pointintervalh",
                          ref_line = NA_real_) {

  # All checks carried out by downstream functions

  s <- summary(x, pars = pars, include = include)
  p <- plot(s, ..., stat = stat, ref_line = ref_line)
  return(p)
}

plot_prior_posterior <- function(x, ...,
                                 prior = NULL,
                                 post_args = list(),
                                 prior_args = list(),
                                 overlay = c("prior", "posterior"),
                                 ref_line = NA_real_) {

  # Checks
  if (!inherits(x, "stan_nma"))
    abort("Not a `stan_nma` object.")

  priors_used <-
    c("intercept"[!is.null(x$priors$prior_intercept)],
      "trt"[!is.null(x$priors$prior_trt)],
      "het"[!is.null(x$priors$prior_het)],
      "reg"[!is.null(x$priors$prior_reg)],
      "aux"[!is.null(x$priors$prior_aux)])

  if (is.null(prior)) {
    prior <- priors_used
  } else if (!rlang::is_character(prior) || !all(prior %in% priors_used)) {
    abort(paste0("`prior` should be a character vector, with elements from ",
                paste(priors_used, collapse = ", ")))
  }

  if (!is.list(post_args))
    abort("`post_args` should be a list of arguments to pass to tidybayes::stat_sample_slabinterval")

  if (!is.list(prior_args))
    abort("`prior_args` should be a list of arguments to pass to tidybayes::stat_dist_slabinterval")

  overlay <- rlang::arg_match(overlay)

  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  # Get prior details
  prior_dat <- vector("list", length(prior))
  for (i in seq_along(prior)) {
    if (prior[i] %in% c("het", "aux")) trunc <- c(0, Inf)
    else trunc <- NULL
    prior_dat[[i]] <- get_tidy_prior(x$priors[[paste0("prior_", prior[i])]], trunc = trunc) %>%
      tibble::add_column(prior = prior[i])
  }

  prior_dat <- dplyr::bind_rows(prior_dat) %>%
    dplyr::mutate(par_base = dplyr::recode(.data$prior,
                                           intercept = "mu",
                                           trt = "d",
                                           het = "tau",
                                           reg = "beta",
                                           aux = "sigma"))


  # Get parameter samples
  pars <- unique(prior_dat$par_base)

  draws <- tibble::as_tibble(as.matrix(x, pars = pars))

  # Transform heterogeneity samples to prior scale (SD, variance, precision)
  if ("het" %in% prior) {
    if (x$priors$prior_het_type == "var") {
      draws$tausq <- draws$tau^2
      draws <- dplyr::select(draws, -.data$tau)
      prior_dat$par_base <- dplyr::recode(prior_dat$par_base, tau = "tausq")
    } else if (x$priors$prior_het_type == "prec") {
      draws$prec <- draws$tau^-2
      draws <- dplyr::select(draws, -.data$tau)
      prior_dat$par_base <- dplyr::recode(prior_dat$par_base, tau = "prec")
    }
  }

  if (packageVersion("tidyr") >= "1.0.0") {
    draws <- tidyr::pivot_longer(draws, cols = dplyr::everything(),
                                 names_to = "parameter", values_to = "value")
  } else {
    draws <- tidyr::gather(key = "parameter",
                           value = "value",
                           dplyr::everything())
  }

  draws$par_base <- stringr::str_remove(draws$parameter, "\\[.*\\]")
  draws$parameter <- forcats::fct_inorder(factor(draws$parameter))

  # Join prior name into posterior
  draws <- dplyr::left_join(draws, prior_dat[, c("par_base", "prior")], by = "par_base")

  # Repeat rows of prior_dat for each corresponding parameter
  prior_dat <- dplyr::left_join(prior_dat,
                                dplyr::distinct(draws, .data$par_base, .data$parameter),
                                by = "par_base")

  # Construct plot

  xlim <- c(min(draws$value), max(draws$value))

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
    ggplot2::coord_cartesian(xlim = xlim)

  g_prior <- rlang::call2(tidybayes::stat_dist_slabinterval,
                          !!! rlang::dots_list(mapping = ggplot2::aes(y = .data$parameter, dist = .data$dist, args = .data$args),
                                               data = prior_dat,
                                               orientation = "horizontal", show_interval = FALSE, normalize = "none",
                                               slab_fill = NA, slab_colour = "black", slab_size = 0.5,
                                               !!! prior_args,
                                               .homonyms = "last"))

  g_post <- rlang::call2(tidybayes::stat_sample_slabinterval,
                         !!! rlang::dots_list(mapping = ggplot2::aes(y = .data$parameter, x = .data$value),
                                              data = draws,
                                              orientation = "horizontal", show_interval = FALSE, normalize = "none",
                                              slab_type = "histogram",
                                              !!! post_args,
                                              .homonyms = "last"))

  if (overlay == "prior") {
    p <- p + eval(g_post) + eval(g_prior)
  } else {
    p <- p + eval(g_prior) + eval(g_post)
  }

  p <- p +
    ggplot2::facet_grid(rows = "prior", scales = "free", space = "free") +
    theme_multinma()

  return(p)
}

#' as.stanfit
#'
#' Attempt to turn an object into a `\link[rstan:stanfit-class]{stanfit}` object.
#'
#' @param x an object
#' @param ... additional arguments
#'
#' @return A `\link[rstan:stanfit-class]{stanfit}` object.
#' @export
#'
#' @examples
as.stanfit <- function(x, ...) {
  UseMethod("as.stanfit")
}

#' @export
#' @rdname as.stanfit
as.stanfit.stan_nma <- function(x, ...) {
  return(x[["stanfit"]])
}

#' @export
#' @rdname as.stanfit
as.stanfit.default <- function(x, ...) {
  abort(glue::glue("Cannot coerce object of class '{class(x)}' to 'stanfit'."))
}

#' as.array
#'
#' Turn a `stan_nma` object into a 3D array \[Iteration, Chain, Parameter\].
#' Enables `\link[bayesplot:bayesplot-package]{bayesplot}` functions to
#' seamlessly work on `stan_nma` objects.
#'
#' @param x an object
#' @param ... additional arguments to [as.array.stanfit]
#'
#' @export
as.array.stan_nma <- function(x, ...) {
  return(as.array(as.stanfit(x), ...))
}

#' as.data.frame
#'
#' Turn a `stan_nma` object into a data frame containing posterior samples (post
#' warm-up) of each parameter.
#'
#' @param x an object
#' @param ... additional arguments to [as.data.frame.stanfit]
#'
#' @export
as.data.frame.stan_nma <- function(x, ...) {
  return(as.data.frame(as.stanfit(x), ...))
}

#' as.matrix
#'
#' Turn a `stan_nma` object into a matrix containing posterior samples (post
#' warm-up) of each parameter.
#'
#' @param x an object
#' @param ... additional arguments to [as.matrix.stanfit]
#'
#' @export
as.matrix.stan_nma <- function(x, ...) {
  return(as.matrix(as.stanfit(x), ...))
}

#' Model comparison using the `loo` package
#'
#' The `\link[loo:loo]{loo()}` and `\link[loo:waic]{waic()}` functions from the `loo`
#' package may be called directly on [stan_nma] and [stan_mlnmr] objects.
#'
#' @param x An object of class [stan_nma] or [stan_mlnmr]
#' @param ... Further arguments to `\link[rstan:stanfit-method-loo]{loo()}` or
#'   `\link[loo:waic]{waic()}`
#'
#' @export
#' @importFrom loo loo
#' @rdname loo
#' @aliases loo
loo.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  return(rstan::loo(sf, ...))
}

#' @export
#' @importFrom loo waic
#' @rdname loo
#' @aliases waic
waic.stan_nma <- function(x, ...) {
  ll <- as.array(x, pars = "log_lik")
  return(loo::waic(ll, ...))
}

#' Matrix of plots for a `stan_nma` object
#'
#' A [pairs()] method for `stan_nma` objects, which calls
#' `\link[rstan:stanfit-method-pairs]{pairs()}` on the underlying `stanfit`
#' object.
#'
#' @param x An object of class `stan_nma`
#' @param ... Other arguments passed to
#'   `\link[rstan:stanfit-method-pairs]{pairs()}`, such as `pars` to select the
#'   parameters to display.
#'
#' @return
#' @export
#'
#' @examples
pairs.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  pairs(sf, ...)
}
