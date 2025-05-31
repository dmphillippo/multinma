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
#'   \code{\link[rstan:stanmodel-method-sampling]{sampling()}} for the model}
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
#'   \item{`basis`}{For `mspline` and `pexp` models, a named list of spline
#'    bases for each study}
#'   }
#'
#' The `stan_mlnmr` sub-class inherits from `stan_nma`, and differs only in the
#' class of the `network` object.
#'
NULL

#' Print `stan_nma` objects
#'
#' @param x A [stan_nma] object
#' @param ... Further arguments passed to \code{\link[rstan:print.stanfit]{print.stanfit()}}
#'
#' @return `x` is returned invisibly.
#'
#' @export
print.stan_nma <- function(x, ...) {
  if (inherits(x$network, "mlnmr_data")) type <- "ML-NMR"
  else type <- "NMA"
  cglue("A {x$trt_effects} effects {type} with a {x$likelihood} likelihood ({x$link} link).")
  if (x$likelihood %in% c("mspline", "pexp")) {
    deg <- switch(x$likelihood,
                  mspline = switch(attr(x$basis[[1]], 'degree'),
                                   "1" = "Piecewise constant",
                                   "2" = "Quadratic M-spline",
                                   "3" = "Cubic M-spline",
                                   paste('Degree', attr(x$basis[[1]], 'degree'), 'M-spline')),
                  pexp = 'Piecewise constant')
    cglue("{deg} baseline hazard with {length(attr(x$basis[[1]], 'knots'))} internal knots.")
  }
  if (x$consistency != "consistency") {
    if (x$consistency == "nodesplit")
      cglue("An inconsistency model ('{x$consistency}') was fitted, splitting the comparison {x$nodesplit[2]} vs. {x$nodesplit[1]}.")
    else
      cglue("An inconsistency model ('{x$consistency}') was fitted.")
  }
  if (!is.null(x$regression)) cglue("Regression model: {rlang::as_label(x$regression)}.")
  if (!is.null(x$aux_regression)) cglue("Auxiliary regression model: {rlang::as_label(x$aux_regression)}.")
  if ((!is.null(x$regression) || !is.null(x$aux_regression)) && !is.null(x$xbar)) {
    cglue("Centred covariates at the following overall mean values:")
    print(x$xbar)
  }
  if (length(setdiff(x$aux_by, ".study"))) cglue("Stratified baseline hazards by {glue::glue_collapse(x$aux_by, sep = ', ', last = ' and ')}.")

  sf <- as.stanfit(x)
  dots <- list(...)
  include <- "pars" %in% names(dots)
  dots <- rlang::dots_list(x = sf,
                           pars = c("log_lik", "resdev",
                                    "fitted_ipd",
                                    "fitted_agd_arm",
                                    "fitted_agd_contrast",
                                    "theta_bar_cum_agd_arm",
                                    "theta_bar_cum_agd_contrast",
                                    "theta2_bar_cum",
                                    "mu", "delta",
                                    if (!is.null(x$aux_regression) &&
                                        length(setdiff(colnames(attr(terms(x$aux_regression), "factor")), ".trt")) > 0) {
                                      if (x$likelihood %in% c("mspline", "pexp")) NULL else "d_aux"
                                    } else {
                                      "beta_aux"
                                    },
                                    "scoef"),
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
#' @param x,object A `stan_nma` object
#' @param ... Additional arguments passed on to other methods
#' @param pars,include See \code{\link[rstan:stanfit-method-extract]{rstan::extract()}}
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
#' ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Summary and plot of all model parameters
#' summary(smk_fit_RE)
#' plot(smk_fit_RE)
#'
#' # Summary and plot of heterogeneity tau only
#' summary(smk_fit_RE, pars = "tau")
#' plot(smk_fit_RE, pars = "tau")
#'
#' # Customising plot output
#' plot(smk_fit_RE,
#'      pars = c("d", "tau"),
#'      stat = "halfeye",
#'      ref_line = 0)
#' }
summary.stan_nma <- function(object, ...,
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
    pars <- c("log_lik", "resdev", "lp__")
    if (has_ipd(object$network)) pars <- c(pars, "fitted_ipd")
    if (has_agd_arm(object$network)) pars <- c(pars, "fitted_agd_arm")
    if (has_agd_contrast(object$network)) pars <- c(pars, "fitted_agd_contrast")
    if (inherits(object, "stan_mlnmr")) {
      if (has_agd_arm(object$network)) pars <- c(pars, "theta_bar_cum_agd_arm")
      if (has_agd_contrast(object$network)) pars <- c(pars, "theta_bar_cum_agd_contrast")
    }
    if (object$likelihood %in% c("bernoulli2", "binomial2")) {
      if (has_agd_arm(object$network)) pars <- c(pars, "theta2_bar_cum")
    }
  } else {
    if (!is.character(pars)) abort("`pars` should be a character vector")
  }

  sims <- as.array(object, pars = pars, include = include)
  sums <- summary_mcmc_array(sims, probs = probs)
  ss <- list(summary = sums, sims = sims)
  class(ss) <- c("nma_parameter_summary", "nma_summary")
  attr(ss, "xlab") <- "Parameter"
  attr(ss, "ylab") <- "Value"
  return(ss)
}

#' @param stat Character string specifying the `ggdist` plot stat to use,
#'   default `"pointinterval"`
#' @param orientation Whether the `ggdist` geom is drawn horizontally
#'   (`"horizontal"`) or vertically (`"vertical"`), default `"horizontal"`
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn
#' @rdname summary.stan_nma
#' @export
plot.stan_nma <- function(x, ...,
                          pars, include,
                          stat = "pointinterval",
                          orientation = c("horizontal", "vertical", "y", "x"),
                          ref_line = NA_real_) {

  # All checks carried out by downstream functions

  s <- summary(x, pars = pars, include = include)
  p <- plot(s, ..., stat = stat, orientation = orientation, ref_line = ref_line)
  return(p)
}

#' Plot prior vs posterior distribution
#'
#' Produce plots comparing the prior and posterior distributions of model
#' parameters.
#'
#' @param x A `stan_nma` object
#' @param ... Additional arguments passed on to methods
#' @param prior Character vector selecting the prior and posterior
#'   distribution(s) to plot. May include `"intercept"`, `"trt"`, `"het"`,
#'   `"reg"`, `"aux"`, `"class_mean"` or `"class_sd"` as appropriate.
#' @param post_args List of arguments passed on to [ggplot2::geom_histogram] to
#'   control plot output for the posterior distribution
#' @param prior_args List of arguments passed on to [ggplot2::geom_path] to
#'   control plot output for the prior distribution. Additionally, `n` controls
#'   the number of points the density curve is evaluated at (default `500`), and
#'   `p_limits` controls the endpoints of the curve as quantiles (default
#'   `c(.001, .999)`).
#' @param overlay String, should prior or posterior be shown on top? Default
#'   `"prior"`.
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn
#'
#' @return A `ggplot` object.
#' @export
#'
#' @details Prior distributions are displayed as lines, posterior distributions
#'   are displayed as histograms.
#'
#' @importFrom truncdist dtrunc ptrunc qtrunc
#'
#' @examples
#' ## Smoking cessation NMA
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Plot prior vs. posterior, by default all parameters are plotted
#' plot_prior_posterior(smk_fit_RE)
#'
#' # Plot prior vs. posterior for heterogeneity SD only
#' plot_prior_posterior(smk_fit_RE, prior = "het")
#'
#' # Customise plot
#' plot_prior_posterior(smk_fit_RE, prior = "het",
#'                      prior_args = list(colour = "darkred", size = 2),
#'                      post_args = list(alpha = 0.6))
#' }
#'
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
      "aux"[!is.null(x$priors$prior_aux)],
      "aux_reg"[!is.null(x$priors$prior_aux_reg)],
      "class_mean"[!is.null(x$priors$prior_class_mean)],
      "class_sd"[!is.null(x$priors$prior_class_sd)])

  if (is.null(prior)) {
    prior <- priors_used
  } else if (!rlang::is_character(prior) || !all(prior %in% priors_used)) {
    abort(paste0("`prior` should be a character vector, with elements from ",
                paste(priors_used, collapse = ", ")))
  }

  if (!is.list(post_args))
    abort("`post_args` should be a list of arguments to pass to ggdist::stat_sample_slabinterval")

  if (!is.list(prior_args))
    abort("`prior_args` should be a list of arguments to pass to ggdist::stat_dist_slabinterval")

  overlay <- rlang::arg_match(overlay)

  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  if (x$likelihood %in% c("mspline", "pexp")) n_scoef <- ncol(x$basis[[1]])

  # Get prior details
  prior_dat <- vector("list", length(prior))
  for (i in seq_along(prior)) {
    if (prior[i] %in% c("het", "aux") || (prior[i] == "aux_reg" && x$likelihood %in% c("mspline", "pexp"))) trunc <- c(0, Inf)
    else trunc <- NULL

    prior_dat[[i]] <- get_tidy_prior(x$priors[[paste0("prior_", prior[i])]], trunc = trunc) %>%
      tibble::add_column(prior = prior[i])

    if (x$likelihood == "gengamma" && prior[i] == "aux") {
      prior_dat[[i]] <-
        dplyr::bind_rows(get_tidy_prior(x$priors$prior_aux$sigma, trunc = trunc),
                         get_tidy_prior(x$priors$prior_aux$k, trunc = trunc)) %>%
        tibble::add_column(prior = c("aux", "aux2"))

    } else {
      prior_dat[[i]] <- get_tidy_prior(x$priors[[paste0("prior_", prior[i])]], trunc = trunc) %>%
        tibble::add_column(prior = prior[i])
    }
  }

  prior_dat <- dplyr::bind_rows(prior_dat) %>%
    dplyr::mutate(par_base = dplyr::recode(.data$prior,
                                           intercept = "mu",
                                           trt = "d",
                                           het = "tau",
                                           reg = "beta",
                                           aux = switch(x$likelihood,
                                                        normal = "sigma",
                                                        ordered = "cc",
                                                        weibull = "shape",
                                                        gompertz = "shape",
                                                        `weibull-aft` = "shape",
                                                        lognormal = "sdlog",
                                                        loglogistic = "shape",
                                                        gamma = "shape",
                                                        gengamma = "sigma",
                                                        mspline = "sigma",
                                                        pexp = "sigma"),
                                           aux2 = switch(x$likelihood,
                                                         gengamma = "k"),
                                           aux_reg = switch(x$likelihood,
                                                            pexp =, mspline = "sigma_beta",
                                                            weibull =, gompertz =, `weibull-aft` =,
                                                            lognormal =, loglogistic =, gamma =,
                                                            gengamma = "beta_aux"),
                                           class_mean = "class_mean",
                                           class_sd = "class_sd"))

  # Add in omega parameter if node-splitting model, which uses prior_trt
  if (inherits(x, "nma_nodesplit")) {
    prior_dat <- dplyr::bind_rows(
      prior_dat,
      dplyr::filter(prior_dat, .data$prior == "trt") %>%
        dplyr::mutate(par_base = "omega")
    )
  }

  # Get parameter samples
  pars <- unique(prior_dat$par_base)

  draws <- tibble::as_tibble(as.matrix(x, pars = pars))

  # Transform heterogeneity samples to prior scale (SD, variance, precision)
  if ("het" %in% prior) {
    if (x$priors$prior_het_type == "var") {
      draws$tausq <- draws$tau^2
      draws <- dplyr::select(draws, -"tau")
      prior_dat$par_base <- dplyr::recode(prior_dat$par_base, tau = "tausq")
    } else if (x$priors$prior_het_type == "prec") {
      draws$prec <- draws$tau^-2
      draws <- dplyr::select(draws, -"tau")
      prior_dat$par_base <- dplyr::recode(prior_dat$par_base, tau = "prec")
    }
  }

  # For ordered likelihood, priors are specified on differences between cutoffs
  if (x$likelihood == "ordered" && "aux" %in% prior) {
    l_cat <- if (has_ipd(x$network)) colnames(x$network$ipd$.r) else colnames(x$network$agd_arm$.r)
    n_cat <- length(l_cat)

    if (n_cat <= 2) {
      draws <- dplyr::select(draws, -dplyr::starts_with("cc["))
      prior_dat <- dplyr::filter(prior_dat, prior != "aux")
    } else {
      for (i in 2:(n_cat-1)) {
        draws <- dplyr::mutate(draws, !! paste0("diff_cc[", l_cat[i+1], " - ", l_cat[i], "]") :=
                                        !! as.symbol(paste0("cc[", l_cat[i+1], "]")) - !! as.symbol(paste0("cc[", l_cat[i], "]")))
      }
      draws <- dplyr::select(draws, -dplyr::starts_with("cc["))
    }

    prior_dat <- dplyr::mutate(prior_dat, par_base = dplyr::recode(.data$par_base, cc = "diff_cc"))
  }

  draws <- tidyr::pivot_longer(draws, cols = dplyr::everything(),
                               names_to = "parameter", values_to = "value")

  draws$par_base <- stringr::str_remove(draws$parameter, "\\[.*\\]")
  draws$parameter <- forcats::fct_inorder(factor(draws$parameter))

  # Join prior name into posterior
  draws <- dplyr::left_join(draws, prior_dat[, c("par_base", "prior")], by = "par_base")

  # Calculate prior density lines
  if (rlang::has_name(prior_args, "p_limits")) {
    p_limits <- prior_args$p_limits
    prior_args <- purrr::list_modify(prior_args, p_limits = purrr::zap())
  } else {
    p_limits <- c(0.001, 0.999)
  }
  if (rlang::has_name(prior_args, "n")) {
    n <- prior_args$n
    prior_args <- purrr::list_modify(prior_args, n = purrr::zap())
  } else {
    n <- 501
  }

  xseq <- dens <- vector("list", nrow(prior_dat))
  for (i in seq_len(nrow(prior_dat))) {
    dist <- prior_dat$dist[[i]]
    args <- prior_dat$args[[i]]

    if (dist == "unif") {
      # lower <- max(args$min, -1e12)
      # upper <- min(args$max, 1e12)

      xseq[[i]] <- c(args$min, args$max)
      if (is.infinite(args$min) || is.infinite(args$max)) {
        dens[[i]] <- c(0, 0)
      } else {
        dens[[i]] <- dunif(xseq[[i]], min = args$min, max = args$max)
      }
    } else {
      lower <- eval(rlang::call2(paste0("q", dist), p = p_limits[1], !!! args))
      upper <- eval(rlang::call2(paste0("q", dist), p = p_limits[2], !!! args))

      xseq[[i]] <- seq(from = lower, to = upper, length.out = n)

      dens[[i]] <- eval(rlang::call2(paste0("d", dist), x = xseq[[i]], !!! args))
    }
  }

  prior_dat <- tibble::add_column(prior_dat, xseq = xseq, dens = dens)
  prior_dat <- tidyr::unnest(prior_dat, c("xseq", "dens"))

  # Repeat rows of prior_dat for each corresponding parameter
  if (packageVersion("dplyr") >= "1.1.1") {
    prior_dat <- dplyr::left_join(prior_dat,
                                  dplyr::distinct(draws, .data$par_base, .data$parameter),
                                  by = "par_base",
                                  relationship = "many-to-many")
  } else {
    prior_dat <- dplyr::left_join(prior_dat,
                                  dplyr::distinct(draws, .data$par_base, .data$parameter),
                                  by = "par_base")
  }

  # Construct plot
  xlim <- c(min(draws$value, 0), max(draws$value))

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
    ggplot2::coord_cartesian(xlim = xlim)

  g_prior <- rlang::call2(ggplot2::geom_line,
                          !!! rlang::dots_list(mapping = ggplot2::aes(x = .data$xseq, y = .data$dens),
                                               data = prior_dat,
                                               !!! prior_args,
                                               .homonyms = "last"))

  g_post <- rlang::call2(ggplot2::geom_histogram,
                         !!! rlang::dots_list(mapping = ggplot2::aes(y = ggplot2::after_stat(.data$density), x = .data$value, group = .data$parameter),
                                              data = draws,
                                              binwidth = function(x) diff(range(x)) / nclass.Sturges(x),
                                              boundary = 0,
                                              position = "identity",
                                              !!! post_args,
                                              .homonyms = "last"))

  if (overlay == "prior") {
    p <- p + eval(g_post) + eval(g_prior)
  } else {
    p <- p + eval(g_prior) + eval(g_post)
  }

  p <- p +
    ggplot2::facet_wrap("parameter", scales = "free") +
    theme_multinma()

  return(p)
}

#' Plot numerical integration error
#'
#' For ML-NMR models, plot the estimated numerical integration error over the
#' entire posterior distribution, as the number of integration points increases.
#' See \insertCite{methods_paper,Phillippo_thesis}{multinma} for details.
#'
#' @param x An object of type `stan_mlnmr`
#' @param ... Additional arguments passed to the `ggdist` plot stat.
#' @param stat Character string specifying the `ggdist` plot stat used to
#'   summarise the integration error over the posterior. Default is `"violin"`,
#'   which is equivalent to `"eye"` with some cosmetic tweaks.
#' @param orientation Whether the `ggdist` geom is drawn horizontally
#'   (`"horizontal"`) or vertically (`"vertical"`), default `"vertical"`
#' @param show_expected_rate Logical, show typical convergence rate \eqn{1/N}?
#'   Default `TRUE`.
#'
#' @details The total number of integration points is set by the `n_int`
#'   argument to [add_integration()], and the intervals at which integration
#'   error is estimated are set by the `int_thin` argument to [nma()]. The
#'   typical convergence rate of Quasi-Monte Carlo integration (as used here) is
#'   \eqn{1/N}, which by default is displayed on the plot output.
#'
#'   The integration error at each thinning interval \eqn{N_\mathrm{thin}}{N_thin} is
#'   estimated for each point in the posterior distribution by subtracting the
#'   final estimate (using all `n_int` points) from the estimate using only the
#'   first \eqn{N_\mathrm{thin}}{N_thin} points.
#'
#' # Note for survival models
#' This function is not supported for survival/time-to-event models. These do
#' not save cumulative integration points for efficiency reasons (both time and
#' memory).
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_network
#' @template ex_plaque_psoriasis_integration
#' @examples \donttest{
#' # Fit the ML-NMR model
#' pso_fit <- nma(pso_net, 
#'                trt_effects = "fixed",
#'                link = "probit",
#'                likelihood = "bernoulli2",
#'                regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
#'                class_interactions = "common",
#'                prior_intercept = normal(scale = 10),
#'                prior_trt = normal(scale = 10),
#'                prior_reg = normal(scale = 10),
#'                init_r = 0.1,
#'                QR = TRUE,
#'                # Set the thinning factor for saving the cumulative results
#'                # (This also sets int_check = FALSE)
#'                int_thin = 8)
#' pso_fit
#'
#' # Plot numerical integration error
#' plot_integration_error(pso_fit)
#' }
plot_integration_error <- function(x, ...,
                                   stat = "violin",
                                   orientation = c("vertical", "horizontal", "x", "y"),
                                   show_expected_rate = TRUE) {
  # Checks
  if (!inherits(x, "stan_mlnmr"))
    abort("Expecting a `stan_mlnmr` object, created by fitting a ML-NMR model with numerical integration using the `nma()` function.")

  if (inherits(x, "stan_nma_surv"))
    abort("Not supported for survival models; cumulative integration points are not saved for efficiency reasons.")

  if (!rlang::is_bool(show_expected_rate))
    abort("`show_expected_rate` must be a logical value, TRUE or FALSE.")

  if (!rlang::is_string(stat))
    abort("`stat` should be a character string specifying the name of a ggdist stat")

  stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

  if (violin <- stat == "violin") {
    stat <- "eye"
  }

  tb_geom <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
                      error = function(err) {
                        abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                                    err, sep = "\n"))
                      })

  orientation <- rlang::arg_match(orientation)
  if (orientation == "x") orientation <- "vertical"
  else if (orientation == "y") orientation <- "horizontal"

  # Is a horizontal geom specified?
  horizontal <- orientation == "horizontal"

  # Get cumulative integration points
  twoparbin <- x$likelihood %in% c("binomial2", "bernoulli2")
  multi <- x$likelihood == "ordered"

  ipars <- c()
  if (has_agd_arm(x$network)) {
    ipars <- c(ipars, "theta_bar_cum_agd_arm")
  }
  if (has_agd_contrast(x$network)) {
    ipars <- c(ipars, "theta_bar_cum_agd_contrast")
  }
  if (twoparbin) {
    ipars <- c(ipars, "theta2_bar_cum")
  }

  if (!all(ipars %in% x$stanfit@sim$pars_oi))
    abort(paste0("Cumulative integration points not saved.\n",
                 "Re-run model with `int_thin > 0` and `int_check = FALSE` to use this feature."))

  int_dat <- as.data.frame(x, pars = ipars) %>%
    dplyr::mutate(.draw = 1:dplyr::n())

  colnames(int_dat) <- stringr::str_remove(colnames(int_dat), "_bar_cum(_agd_arm|_agd_contrast)?")

  n_int <- x$network$n_int

  rx <- if (multi) "^(theta2?)\\[(.+): (.+), ([0-9]+), (.+)\\]$" else "^(theta2?)\\[(.+): (.+), ([0-9]+)\\]$"

  int_dat <- tidyr::pivot_longer(int_dat, cols = -dplyr::one_of(".draw"),
                                 names_pattern = rx,
                                 names_to = if (multi) c("parameter", "study", "treatment", "n_int", "category") else c("parameter", "study", "treatment", "n_int"),
                                 names_transform = list(n_int = as.integer),
                                 values_to = "value")

  int_dat$study <- factor(int_dat$study, levels = levels(x$network$studies))
  int_dat$treatment <- factor(int_dat$treatment, levels = levels(x$network$treatments))
  if (multi) int_dat$category <- forcats::fct_inorder(factor(int_dat$category))

  # Estimate integration error by subtracting final value
  int_dat <- dplyr::left_join(dplyr::filter(int_dat, .data$n_int != max(.data$n_int)),
                              dplyr::filter(int_dat, .data$n_int == max(.data$n_int)) %>%
                                dplyr::rename(final_value = "value") %>%
                                dplyr::select(-"n_int"),
                              by = if (multi)  c("parameter", "study", "treatment", "category", ".draw") else c("parameter", "study", "treatment", ".draw")) %>%
    dplyr::mutate(diff = .data$value - .data$final_value)

  int_thin <- min(int_dat$n_int)

  # Reference convergence rates
  if (show_expected_rate) {
    conv_dat <- dplyr::bind_rows(
      tibble::tibble(n_int = seq(1, n_int, length.out = 501), diff = n_int^-1, group = "pos"),
      tibble::tibble(n_int = seq(1, n_int, length.out = 501), diff = -n_int^-1, group = "neg"))
  }

  # Create plot
  if (horizontal) {
    p <- ggplot2::ggplot(int_dat, ggplot2::aes(y = .data$n_int, x = .data$diff)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      ggplot2::xlab("Estimated integration error") +
      ggplot2::ylab("Number of integration points") +
      ggplot2::coord_cartesian(xlim = range(int_dat$diff))
  } else {
    p <- ggplot2::ggplot(int_dat, ggplot2::aes(x = .data$n_int, y = .data$diff)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::ylab("Estimated integration error") +
      ggplot2::xlab("Number of integration points") +
      ggplot2::coord_cartesian(ylim = range(int_dat$diff))
  }

  if (show_expected_rate) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(group = .data$group),
                         data = conv_dat,
                         colour = "grey60", linetype = 2)
  }

  # Custom format for violin plot
  v_args <-
    if (violin) {
      if (twoparbin) {
        list(point_interval = NULL,
             alpha = 0.8,
             slab_size = 0.5)
        } else {
          list(point_interval = NULL,
               slab_colour = "black",
               slab_size = 0.5)
        }
    } else list()

  if (twoparbin) {
    p <- p +
      do.call(tb_geom,
              args = rlang::dots_list(ggplot2::aes(colour = .data$parameter,
                                                   fill = .data$parameter,
                                                   slab_colour = .data$parameter),
                                      orientation = orientation,
                                      ...,
                                      !!! v_args,
                                      position = if (!horizontal) ggplot2::position_dodge(width = int_thin / 10) else "identity",
                                      .homonyms = "first")) +
      ggplot2::scale_colour_manual("Parameter",
                                   values = c(theta = "#113259", theta2 = "#55A480"),
                                   labels = function(x) dplyr::recode(x,
                                     theta = expression(bar(p)),
                                     theta2 = expression(bar(p)^2)),
                                   aesthetics = c("colour", "fill", "slab_colour"))
  } else {
    p <- p +
      do.call(tb_geom,
              args = rlang::dots_list(orientation = orientation, ..., !!! v_args, .homonyms = "first"))
  }

  if (multi) {
    p <- p + ggplot2::facet_grid(study + treatment ~ category)
  } else {
    p <- p + ggplot2::facet_wrap(~ study + treatment)
  }

  p <- p + theme_multinma()

  return(p)
}

#' as.stanfit
#'
#' Attempt to turn an object into a \code{\link[rstan:stanfit-class]{stanfit}} object.
#'
#' @param x an object
#' @param ... additional arguments
#'
#' @return A \code{\link[rstan:stanfit-class]{stanfit}} object.
#' @export
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

#' Convert samples into arrays, matrices, or data frames
#'
#' Samples (post warm-up) from a `stan_nma` model object can be coerced into an
#' array, matrix, or data frame.
#'
#' @param x A `stan_nma` object
#' @param ... Additional arguments passed to \code{\link[rstan:as.array.stanfit]{as.array.stanfit()}}
#' @param pars Optional character vector of parameter names to include in output. If not specified, all parameters are used.
#' @param include Logical, are parameters in `pars` to be included (`TRUE`, default) or excluded (`FALSE`)?
#'
#' @return The `as.array()` method produces a 3D array \[Iteration, Chain,
#'   Parameter\] containing posterior samples of each parameter (as class
#'   [mcmc_array]). This has the side effect of enabling
#'   \code{\link[bayesplot:bayesplot-package]{bayesplot}} functions to
#'   seamlessly work on `stan_nma` objects.
#'
#'   The `as.data.frame()` method produces a data frame containing posterior
#'   samples of each parameter, combined over all chains.
#'
#'   The `as.matrix()` method produces a matrix containing posterior samples of
#'   each parameter, combined over all chains.
#'
#' @export
as.array.stan_nma <- function(x, ..., pars, include = TRUE) {
  if (!rlang::is_bool(include))
    abort("`include` must be TRUE or FALSE.")
  if (!missing(pars)) {
    if (!is.character(pars))
      abort("`pars` must be a character vector of parameters to include (or missing).")

    allpars <- c(x$stanfit@sim$pars_oi, x$stanfit@sim$fnames_oi)
    badpars <- setdiff(pars, allpars)
    if (length(badpars))
      abort(glue::glue("No parameter{if (length(badpars) > 1) 's' else ''} ",
                       glue::glue_collapse(glue::double_quote(badpars), sep = ", ", last = " or "), "."))

    # Extract from stanfit only parameters represented in pars
    if (include) {
      par_base <- stringr::str_remove(pars, "\\[.*$")
    } else {
      par_base <- setdiff(x$stanfit@sim$pars_oi, pars)
    }
    a <- as.array(as.stanfit(x), pars = par_base, ...)

    # Get parameter indices, respecting order of `pars`
    par_regex <- paste0("^\\Q", pars, "\\E(\\[|$)")
    par_select <- unique(unlist(lapply(par_regex, grep, x = dimnames(a)[[3]],
                                       perl = TRUE, invert = !include)))
    out <- a[ , , par_select, drop = FALSE]

  } else {
    out <- as.array(as.stanfit(x), ...)
  }

  class(out) <- c("mcmc_array", "array")

  return(out)
}

#' @rdname as.array.stan_nma
#' @export
as.data.frame.stan_nma <- function(x, ..., pars, include = TRUE) {
  return(as.data.frame(as.matrix(x, ..., pars = pars, include = include)))
}

#' @rdname as.array.stan_nma
#' @export
as_tibble.stan_nma <- function(x, ..., pars, include = TRUE) {
  return(tibble::as_tibble(as.matrix(x, ..., pars = pars, include = include)))
}

#' @rdname as.array.stan_nma
#' @method as.tibble stan_nma
#' @export
as.tibble.stan_nma <- function(x, ..., pars, include = TRUE) {
  return(tibble::as_tibble(as.matrix(x, ..., pars = pars, include = include)))
}

#' @rdname as.array.stan_nma
#' @export
as.matrix.stan_nma <- function(x, ..., pars, include = TRUE) {
  a <- as.array(x, ..., pars = pars, include = include)
  names_a <- dimnames(a)
  dim_a <- dim(a)
  dim(a) <- c(dim_a[1] * dim_a[2], dim_a[3])
  dimnames(a) <- names_a[-2]
  class(a) <- "matrix"
  return(a)
}

#' Model comparison using the `loo` package
#'
#' The \code{\link[loo:loo]{loo()}} and \code{\link[loo:waic]{waic()}} functions from the `loo`
#' package may be called directly on [stan_nma] and [stan_mlnmr] objects.
#'
#' @param x An object of class [stan_nma] or [stan_mlnmr]
#' @param ... Further arguments to \code{\link[rstan:stanfit-method-loo]{loo()}} or
#'   \code{\link[loo:waic]{waic()}}
#'
#' @rdname loo
#' @aliases loo
#' @method loo stan_nma
# Dynamically exported, see zzz.R
loo.stan_nma <- function(x, ...) {
  sf <- as.stanfit(x)
  return(rstan::loo(sf, ...))
}


#' @rdname loo
#' @aliases waic
#' @method waic stan_nma
# Dynamically exported, see zzz.R
waic.stan_nma <- function(x, ...) {
  ll <- as.array(x, pars = "log_lik")
  return(loo::waic(ll, ...))
}

#' Matrix of plots for a `stan_nma` object
#'
#' A [pairs()] method for `stan_nma` objects, which calls
#' \code{\link[bayesplot:MCMC-scatterplots]{bayesplot::mcmc_pairs()}} on the
#' underlying `stanfit` object.
#'
#' @param x An object of class `stan_nma`
#' @param ... Other arguments passed to
#'   \code{\link[bayesplot:MCMC-scatterplots]{bayesplot::mcmc_pairs()}}
#' @param pars Optional character vector of parameter names to include in
#'   output. If not specified, all parameters are used.
#' @param include Logical, are parameters in `pars` to be included (`TRUE`,
#'   default) or excluded (`FALSE`)?
#'
#' @return A grid of ggplot objects produced by
#'   \code{\link[bayesplot:MCMC-scatterplots]{bayesplot::mcmc_pairs()}}.
#' @export
#'
#' @examples \dontrun{
#' ## Parkinson's mean off time reduction
#' park_net <- set_agd_arm(parkinsons,
#'                         study = studyn,
#'                         trt = trtn,
#'                         y = y,
#'                         se = se,
#'                         sample_size = n)
#'
#' # Fitting a RE model
#' park_fit_RE <- nma(park_net,
#'                    trt_effects = "random",
#'                    prior_intercept = normal(scale = 100),
#'                    prior_trt = normal(scale = 100),
#'                    prior_het = half_normal(scale = 5))
#'
#' # We see a small number of divergent transition errors
#' # These do not go away entirely when adapt_delta is increased
#'
#' # Try to diagnose with a pairs plot
#' pairs(park_fit_RE, pars = c("mu[4]", "d[3]", "delta[4: 3]", "tau"))
#'
#' # Transforming tau onto log scale
#' pairs(park_fit_RE, pars = c("mu[4]", "d[3]", "delta[4: 3]", "tau"),
#'       transformations = list(tau = "log"))
#'
#' # The divergent transitions occur in the upper tail of the heterogeneity
#' # standard deviation. In this case, with only a small number of studies, there
#' # is not very much information to estimate the heterogeneity standard deviation
#' # and the prior distribution may be too heavy-tailed. We could consider a more
#' # informative prior distribution for the heterogeneity variance to aid
#' # estimation.
#' }
#'
pairs.stan_nma <- function(x, ..., pars, include = TRUE) {
  sf <- as.stanfit(x)
  post_array <- as.array(x, pars = pars, include = include)

  max_td <- sf@stan_args[[1]]$control$max_treedepth
  if (is.null(max_td)) max_td <- 10

  args <- rlang::dots_list(x = post_array,
                           np = bayesplot::nuts_params(sf),
                           lp = bayesplot::log_posterior(sf),
                           max_treedepth = max_td,
                           ...,
                           condition = bayesplot::pairs_condition(nuts = "accept_stat__"),
                           .homonyms = "first")

  thm <- ggplot2::theme_set(theme_multinma())
  out <- do.call(bayesplot::mcmc_pairs, args = args)
  ggplot2::theme_set(thm)
  return(out)
}
