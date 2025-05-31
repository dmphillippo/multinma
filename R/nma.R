#' Network meta-analysis models
#'
#' The `nma` function fits network meta-analysis and (multilevel) network
#' meta-regression models in Stan.
#'
#' @param network An `nma_data` object, as created by the functions `set_*()`,
#'   `combine_network()`, or `add_integration()`
#' @param consistency Character string specifying the type of (in)consistency
#'   model to fit, either `"consistency"`, `"ume"`, or `"nodesplit"`
#' @param trt_effects Character string specifying either `"fixed"` or `"random"`
#'   effects
#' @param regression A one-sided model formula, specifying the prognostic and
#'   effect-modifying terms for a regression model. Any references to treatment
#'   should use the `.trt` special variable, for example specifying effect
#'   modifier interactions as `variable:.trt` (see details).
#' @param class_interactions Character string specifying whether effect modifier
#'   interactions are specified as `"common"`, `"exchangeable"`, or
#'   `"independent"`.
#' @param class_effects Character string specifying a model for treatment class effects,
#'   either `"independent"` (the default), `"exchangeable"`, or `"common"`.
#' @param class_sd Character string specifying whether the class standard deviations in a
#'   class effects model should be `"independent"` (i.e. separate for each class, the default),
#'   or `"common"` (i.e. shared across all classes). Alternatively this can be a list of
#'   character vectors, each of which describe a set classes for which to share a common class SD;
#'   any list names will be used to name the output parameters, otherwise the name will be taken
#'   from the first class in each set.
#' @param likelihood Character string specifying a likelihood, if unspecified
#'   will be inferred from the data (see details)
#' @param link Character string specifying a link function, if unspecified will
#'   default to the canonical link (see details)
#' @param ... Further arguments passed to
#'   \code{\link[rstan:stanmodel-method-sampling]{sampling()}}, such as `iter`,
#'   `chains`, `cores`, etc.
#' @param nodesplit For `consistency = "nodesplit"`, the comparison(s) to split
#'   in the node-splitting model(s). Either a length 2 vector giving the
#'   treatments in a single comparison, or a 2 column data frame listing
#'   multiple treatment comparisons to split in turn. By default, all possible
#'   comparisons will be chosen (see [get_nodesplits()]).
#' @param prior_intercept Specification of prior distribution for the intercept
#' @param prior_trt Specification of prior distribution for the treatment
#'   effects
#' @param prior_het Specification of prior distribution for the heterogeneity
#'   (if `trt_effects = "random"`)
#' @param prior_het_type Character string specifying whether the prior
#'   distribution `prior_het` is placed on the heterogeneity standard deviation
#'   \eqn{\tau} (`"sd"`, the default), variance \eqn{\tau^2} (`"var"`), or
#'   precision \eqn{1/\tau^2} (`"prec"`).
#' @param prior_reg Specification of prior distribution for the regression
#'   coefficients (if `regression` formula specified)
#' @param prior_aux Specification of prior distribution for the auxiliary
#'   parameter, if applicable (see details). For `likelihood = "gengamma"` this
#'   should be a list of prior distributions with elements `sigma` and `k`.
#' @param prior_aux_reg Specification of prior distribution for the auxiliary
#'   regression parameters, if `aux_regression` is specified (see details).
#' @param aux_by Vector of variable names listing the variables to stratify the
#'   auxiliary parameters by. Currently only used for survival models, see
#'   details. Cannot be used with `aux_regression`.
#' @param aux_regression A one-sided model formula giving a regression model for
#'   the auxiliary parameters. Currently only used for survival models, see
#'   details. Cannot be used with `aux_by`.
#' @param prior_class_mean Specification of prior distribution for the
#'   treatment class means (if `class_effects = "exchangeable"`).
#' @param prior_class_sd Specification of prior distribution for the
#'   treatment class standard deviations (if `class_effects = "exchangeable"`).
#' @param QR Logical scalar (default `FALSE`), whether to apply a QR
#'   decomposition to the model design matrix
#' @param center Logical scalar (default `TRUE`), whether to center the
#'   (numeric) regression terms about the overall means
#' @param adapt_delta See [adapt_delta] for details
#' @param int_thin A single integer value, the thinning factor for returning
#'   cumulative estimates of integration error. Saving cumulative estimates is
#'   disabled by `int_thin = 0`, which is the default.
#' @param int_check Logical, check sufficient accuracy of numerical integration
#'   by fitting half of the chains with `n_int/2`? When `TRUE`, `Rhat` and
#'   `n_eff` diagnostic warnings will be given if numerical integration has not
#'   sufficiently converged (suggesting increasing `n_int` in
#'   [add_integration()]). Default `TRUE`, but disabled (`FALSE`) when
#'   `int_thin > 0`.
#' @param mspline_degree Non-negative integer giving the degree of the M-spline
#'   polynomial for `likelihood = "mspline"`. Piecewise exponential hazards
#'   (`likelihood = "pexp"`) are a special case with `mspline_degree = 0`.
#' @param n_knots For `mspline` and `pexp` likelihoods, a non-negative integer
#'   giving the number of internal knots for partitioning the baseline hazard
#'   into intervals. The knot locations within each study will be determined by
#'   the corresponding quantiles of the observed event times, plus boundary
#'   knots at the earliest entry time (0 with no delayed entry) and the maximum
#'   event/censoring time. For example, with `n_knots = 3`, the internal knot
#'   locations will be at the 25%, 50%, and 75% quantiles of the observed event
#'   times. The default is `n_knots = 7`; overfitting is avoided by shrinking
#'   towards a constant hazard with a random walk prior (see details). If
#'   `aux_regression` is specified then a single set of knot locations will be
#'   calculated across all studies in the network. See [make_knots()] for more
#'   details on the knot positioning algorithms. Ignored when `knots` is
#'   specified.
#' @param knots For `mspline` and `pexp` likelihoods, a named list of numeric
#'   vectors of knot locations (including boundary knots) for each of the
#'   studies in the network. Currently, each vector must have the same length
#'   (i.e. each study must use the same number of knots). Alternatively, a
#'   single numeric vector of knot locations can be provided which will be
#'   shared across all studies in the network. If unspecified (the default), the
#'   knots will be chosen based on `n_knots` as described above. If
#'   `aux_regression` is specified then `knots` should be a single numeric
#'   vector of knot locations which will be shared across all studies in the
#'   network. [make_knots()] can be used to help specify `knots` directly, or to
#'   investigate knot placement prior to model fitting.
#' @param mspline_basis Instead of specifying `mspline_degree` and `n_knots` or
#'   `knots`, a named list of M-spline bases (one for each study) can be
#'   provided with `mspline_basis` which will be used directly. In this case,
#'   all other M-spline options will be ignored.
#'
#' @details When specifying a model formula in the `regression` argument, the
#'   usual formula syntax is available (as interpreted by [model.matrix()]). The
#'   only additional requirement here is that the special variable `.trt` should
#'   be used to refer to treatment. For example, effect modifier interactions
#'   should be specified as `variable:.trt`. Prognostic (main) effects and
#'   interactions can be included together compactly as `variable*.trt`, which
#'   expands to `variable + variable:.trt` (plus `.trt`, which is already in the
#'   NMA model).
#'
#'   For the advanced user, the additional specials `.study` and `.trtclass` are
#'   also available, and refer to studies and (if specified) treatment classes
#'   respectively. When node-splitting models are fitted (`consistency =
#'   "nodesplit"`) the special `.omega` is available, indicating the arms to
#'   which the node-splitting inconsistency factor is added.
#'
#'   See \code{\link[multinma:priors]{?priors}} for details on prior
#'   specification. Default prior distributions are available, but may not be
#'   appropriate for the particular setting and will raise a warning if used. No
#'   attempt is made to tailor these defaults to the data provided. Please
#'   consider appropriate prior distributions for the particular setting,
#'   accounting for the scales of outcomes and covariates, etc. The function
#'   [plot_prior_posterior()] may be useful in examining the influence of the
#'   chosen prior distributions on the posterior distributions, and the
#'   \code{\link[multinma:summary.nma_prior]{summary()}} method for `nma_prior`
#'   objects prints prior intervals.
#'
#' @section Likelihoods and link functions:
#'   Currently, the following likelihoods and link functions are supported for
#'   each data type:
#'
#'   | \strong{Data type} | \strong{Likelihood}   | \strong{Link function} |
#'   |--------------------|-----------------------|------------------------|
#'   | \strong{Binary}    | `bernoulli`, `bernoulli2`| `logit`, `probit`, `cloglog`
#'   | \strong{Count}     | `binomial`, `binomial2`  | `logit`, `probit`, `cloglog`
#'   | \strong{Rate}      | `poisson`    | `log`
#'   | \strong{Continuous}| `normal`     | `identity`, `log`
#'   | \strong{Ordered}   | `ordered`    | `logit`, `probit`, `cloglog`
#'   | \strong{Survival}  | `exponential`, `weibull`, `gompertz`, `exponential-aft`, `weibull-aft`, `lognormal`, `loglogistic`, `gamma`, `gengamma`, `mspline`, `pexp` | `log`
#'
#'   The `bernoulli2` and `binomial2` likelihoods correspond to a two-parameter
#'   Binomial likelihood for arm-based AgD, which more closely matches the
#'   underlying Poisson Binomial distribution for the summarised aggregate
#'   outcomes in a ML-NMR model than the typical (one parameter) Binomial
#'   distribution \insertCite{@see @methods_paper}{multinma}.
#'
#'   When a `cloglog` link is used, including an offset for log follow-up time
#'   (i.e. `regression = ~offset(log(time))`) results in a model on the log
#'   hazard \insertCite{@see @TSD2}{multinma}.
#'
#'   For survival data, all accelerated failure time models (`exponential-aft`,
#'   `weibull-aft`, `lognormal`, `loglogistic`, `gamma`, `gengamma`) are
#'   parameterised so that the treatment effects and any regression parameters
#'   are log Survival Time Ratios (i.e. a coefficient of \eqn{\log(2)} means
#'   that the treatment or covariate is associated with a doubling of expected
#'   survival time). These can be converted to log Acceleration Factors using
#'   the relation \eqn{\log(\mathrm{AF}) = -\log(\mathrm{STR})} (or equivalently
#'   \eqn{\mathrm{AF} = 1/\mathrm{STR}}).
#'
#'   Further details on each likelihood and link function are given by
#'   \insertCite{TSD2;textual}{multinma}.
#'
#'
#' @section Auxiliary parameters:
#'   Auxiliary parameters are only present in the following models.
#'
#'   ## Normal likelihood with IPD
#'   When a Normal likelihood is fitted to IPD, the auxiliary parameters are the
#'   arm-level standard deviations \eqn{\sigma_{jk}} on treatment \eqn{k} in
#'   study \eqn{j}.
#'
#'   ## Ordered multinomial likelihood
#'   When fitting a model to \eqn{M} ordered outcomes, the auxiliary parameters
#'   are the latent cutoffs between each category, \eqn{c_0 < c_1 < \dots <
#'   c_M}. Only \eqn{c_2} to \eqn{c_{M-1}} are estimated; we fix \eqn{c_0 =
#'   -\infty}, \eqn{c_1 = 0}, and \eqn{c_M = \infty}. When specifying priors for
#'   these latent cutoffs, we choose to specify priors on the *differences*
#'   \eqn{c_{m+1} - c_m}. Stan automatically truncates any priors so that the
#'   ordering constraints are satisfied.
#'
#'   ## Survival (time-to-event) likelihoods
#'   All survival likelihoods except the `exponential` and `exponential-aft`
#'   likelihoods have auxiliary parameters. These are typically study-specific
#'   shape parameters \eqn{\gamma_j>0}, except for the `lognormal` likelihood
#'   where the auxiliary parameters are study-specific standard deviations on
#'   the log scale \eqn{\sigma_j>0}.
#'
#'   The `gengamma` likelihood has two sets of auxiliary parameters,
#'   study-specific scale parameters \eqn{\sigma_j>0} and shape parameters
#'   \eqn{k_j}, following the parameterisation of
#'   \insertCite{Lawless1980;textual}{multinma}, which permits a range of
#'   behaviours for the baseline hazard including increasing, decreasing,
#'   bathtub and arc-shaped hazards. This parameterisation is related to that
#'   discussed by \insertCite{Cox2007;textual}{multinma} and implemented in the
#'   `flexsurv` package with \eqn{Q = k^{-0.5}}. The parameterisation used here
#'   effectively bounds the shape parameter \eqn{k} away from numerical
#'   instabilities as \eqn{k \rightarrow \infty} (i.e. away from \eqn{Q
#'   \rightarrow 0}, the log-Normal distribution) via the prior distribution.
#'   Implicitly, this parameterisation is restricted to \eqn{Q > 0} and so
#'   certain survival distributions like the inverse-Gamma and inverse-Weibull
#'   are not part of the parameter space; however, \eqn{Q > 0} still encompasses
#'   the other survival distributions implemented in this package.
#'
#'   For the `mspline` and `pexp` likelihoods, the auxiliary parameters are the
#'   spline coefficients for each study. These form a unit simplex (i.e. lie
#'   between 0 and 1, and sum to 1), and are given a random walk prior
#'   distribution. `prior_aux` specifies the hyperprior on the random walk
#'   standard deviation \eqn{\sigma} which controls the level of smoothing of
#'   the baseline hazard, with \eqn{\sigma = 0} corresponding to a constant
#'   baseline hazard.
#'
#'   The auxiliary parameters can be stratified by additional factors through
#'   the `aux_by` argument. For example, to allow the shape of the baseline
#'   hazard to vary between treatment arms as well as studies, use `aux_by =
#'   c(".study", ".trt")`. (Technically, `.study` is always included in the
#'   stratification even if omitted from `aux_by`, but we choose here to make
#'   the stratification explicit.) This is a common way of relaxing the
#'   proportional hazards assumption. The default is equivalent to `aux_by =
#'   ".study"` which stratifies the auxiliary parameters by study, as described
#'   above.
#'
#'   A regression model may be specified on the auxiliary parameters using
#'   `aux_regression`. This is useful if we wish to model departures from
#'   non-proportionality, rather than allowing the baseline hazards to be
#'   completely independent using `aux_by`. This is necessary if absolute
#'   predictions (e.g. survival curves) are required in a population for
#'   unobserved combinations of covariates; for example, if `aux_by = .trt` then
#'   absolute predictions may only be produced for the observed treatment arms
#'   in each study population, whereas if `aux_regression = ~.trt` then absolute
#'   predictions can be produced for all treatments in any population. For
#'   `mspline` and `pexp` likelihoods, the regression coefficients are smoothed
#'   over time using a random walk prior to avoid overfitting: `prior_aux_reg`
#'   specifies the hyperprior for the random walk standard deviation. For other
#'   parametric likelihoods, `prior_aux_reg` specifies the prior for the
#'   auxiliary regression coefficients.
#'
#' @return `nma()` returns a [stan_nma] object, except when `consistency =
#'   "nodesplit"` when a [nma_nodesplit] or [nma_nodesplit_df] object is
#'   returned. `nma.fit()` returns a \code{\link[rstan:stanfit-class]{stanfit}}
#'   object.
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' ## Smoking cessation NMA
#' @template ex_smoking_network
#' @template ex_smoking_nma_fe
#' @template ex_smoking_nma_re
#' @template ex_smoking_nma_re_ume
#' @template ex_smoking_nma_re_nodesplit
#' @examples \donttest{
#' # Summarise the node-splitting results
#' summary(smk_fit_RE_nodesplit)
#' }
#'
#' ## Plaque psoriasis ML-NMR
#' @template ex_plaque_psoriasis_network
#' @template ex_plaque_psoriasis_integration
#' @template ex_plaque_psoriasis_mlnmr
#'
#' @examples
#' ## Newly-diagnosed multiple myeloma NMA
#' @template ex_ndmm_network
#' @template ex_ndmm
nma <- function(network,
                consistency = c("consistency", "ume", "nodesplit"),
                trt_effects = c("fixed", "random"),
                regression = NULL,
                class_interactions = c("common", "exchangeable", "independent"),
                class_effects = c("independent", "common", "exchangeable"),
                class_sd =  c("independent", "common"),
                likelihood = NULL,
                link = NULL,
                ...,
                nodesplit = get_nodesplits(network, include_consistency = TRUE),
                prior_intercept = .default(normal(scale = 100)),
                prior_trt = .default(normal(scale = 10)),
                prior_het = .default(half_normal(scale = 5)),
                prior_het_type = c("sd", "var", "prec"),
                prior_reg = .default(normal(scale = 10)),
                prior_aux = .default(),
                prior_aux_reg = .default(),
                prior_class_mean = .default(normal(scale = 10)),
                prior_class_sd = .default(half_normal(scale = 5)),
                aux_by = NULL,
                aux_regression = NULL,
                QR = FALSE,
                center = TRUE,
                adapt_delta = NULL,
                int_thin = 0,
                int_check = TRUE,
                mspline_degree = 3,
                n_knots = 7,
                knots = NULL,
                mspline_basis = NULL) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # Check model arguments
  consistency <- rlang::arg_match(consistency)
  if (length(consistency) > 1) abort("`consistency` must be a single string.")
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")
  class_effects <- rlang::arg_match(class_effects)
  if (length(class_effects) > 1) abort("`class_effects` must be a single string.")

  # Check class_effects and network classes
  if (class_effects != "independent") {
    if (is.null(network$classes)) {
      abort(paste("Setting `class_effects` requires treatment classes to be specified in the network.",
                  "See set_*() argument `trt_class`.", sep = "\n"))
    }
  }



  if (class_effects == "common") {
    # Overwrite treatments with class variables
    if (has_ipd(network)) {
      network$ipd$.trt <- network$ipd$.trtclass
    }
    if (has_agd_arm(network)) {
      network$agd_arm$.trt <- network$agd_arm$.trtclass
    }
    if (has_agd_contrast(network)) {
      network$agd_contrast$.trt <- network$agd_contrast$.trtclass
    }

    # Set the network treatments vector
    network$treatments <- if (.is_default(network$treatments)) {
        .default(factor(levels(network$classes), levels = levels(network$classes)))
      } else {
        factor(levels(network$classes), levels = levels(network$classes))
      }

    # Set network classes vector
    network$classes <- network$treatments
  }

  # Check class_sd
  if (is.list(class_sd)) {
    # Check that all classes listed in 'class_sd' are in 'network$classes'
    if (!all(unlist(class_sd) %in% network$classes)) {
      stop("Some classes listed in 'class_sd' are not present in the network.")
    }


    # Check that all the collapsed classes are distinct and don't share a class
    flattened_classes <- unlist(class_sd)
    if (length(flattened_classes) != length(unique(flattened_classes))) {
      stop("Some classes are listed in more than one shared standard deviation group in 'class_sd'")
    }
  } else {
    class_sd <- rlang::arg_match(class_sd)
    if (length(class_sd) > 1) abort("`class_sd` must be a single string.")
  }

  if (consistency == "nodesplit") {

    lvls_trt <- levels(network$treatments)
    nodesplit_include_consistency <- FALSE

    if (is.data.frame(nodesplit)) { # Data frame listing comparisons to split
      if (ncol(nodesplit) != 2)
        abort("The data frame passed to `nodesplit` should have two columns.")

      nodesplit <- tibble::as_tibble(nodesplit)
      colnames(nodesplit) <- c("trt1", "trt2")

      # NA rows indicate include_consistency = TRUE, filter these out
      if (any(is.na(nodesplit[,1]) & is.na(nodesplit[,2]))) {
        nodesplit_include_consistency <- TRUE
        nodesplit <- dplyr::filter(nodesplit, !is.na(.data$trt1) & !is.na(.data$trt2))
      }

      if (nrow(nodesplit) == 0) {
        abort("No comparisons to node-split.")
      }

      nodesplit$trt1 <- as.character(nodesplit$trt1)
      nodesplit$trt2 <- as.character(nodesplit$trt2)

      if (!all(unlist(nodesplit) %in% lvls_trt))
        abort(sprintf("All comparisons in `nodesplit` should match two treatments in the network.\nSuitable values are: %s",
                      ifelse(length(lvls_trt) <= 5,
                             paste0(lvls_trt, collapse = ", "),
                             paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

      if (any(nodesplit[,1] == nodesplit[,2]))
        abort("`nodesplit` comparison cannot be the same treatment against itself.")

      # Check valid nodesplit - must have both direct and indirect evidence
      ns_check <- dplyr::rowwise(nodesplit) %>%
        dplyr::mutate(direct = has_direct(network, .data$trt1, .data$trt2),
                      indirect = has_indirect(network, .data$trt1, .data$trt2),
                      valid = .data$direct && .data$indirect)

      if (any(!ns_check$valid)) {
        ns_valid <- dplyr::filter(ns_check, .data$valid) %>%
          dplyr::ungroup() %>%
          dplyr::select("trt1", "trt2")

        ns_invalid <- dplyr::filter(ns_check, !.data$valid) %>%
          dplyr::mutate(comparison = paste(.data$trt1, .data$trt2, sep = " vs. "))

        if (nrow(ns_valid)) {
          warn(glue::glue(
            "Ignoring node-split comparisons without both both direct and independent indirect evidence: ",
            glue::glue_collapse(ns_invalid$comparison, sep = ", ", width = 100), "."
            ))

          nodesplit <- ns_valid
        } else {
          abort("No valid comparisons for node-splitting given in `nodesplit`.\n Comparisons must have both direct and independent indirect evidence for node-splitting.")
        }

      }

      # Store comparisons as factors, in increasing order (i.e. trt1 < trt2)
      nodesplit$trt1 <- factor(nodesplit$trt1, levels = lvls_trt)
      nodesplit$trt2 <- factor(nodesplit$trt2, levels = lvls_trt)

      for (i in 1:nrow(nodesplit)) {
        if (as.numeric(nodesplit$trt1[i]) > as.numeric(nodesplit$trt2[i])) {
          nodesplit[i, ] <- rev(nodesplit[i, ])
        }
      }

      # Iteratively call node-splitting models
      n_ns <- nrow(nodesplit) + nodesplit_include_consistency
      ns_fits <- vector("list", n_ns)

      ns_arglist <- list(network = network,
                         consistency = "nodesplit",
                         trt_effects = trt_effects,
                         class_effects = class_effects,
                         class_sd = class_sd,
                         regression = regression,
                         likelihood = likelihood,
                         link = link,
                         ...,
                         prior_intercept = prior_intercept,
                         prior_trt = prior_trt,
                         prior_het = prior_het,
                         prior_het_type = prior_het_type,
                         prior_reg = prior_reg,
                         prior_aux = prior_aux,
                         prior_aux_reg = prior_aux_reg,
                         aux_by = aux_by,
                         aux_regression = aux_regression,
                         prior_class_mean = prior_class_mean,
                         prior_class_sd = prior_class_sd,
                         QR = QR,
                         center = center,
                         adapt_delta = adapt_delta,
                         int_thin = int_thin,
                         mspline_degree = mspline_degree,
                         n_knots = n_knots,
                         knots = knots,
                         mspline_basis = mspline_basis)

      if (!missing(class_interactions)) ns_arglist$class_interactions <- class_interactions


      for (i in 1:nrow(nodesplit)) {

        inform(glue::glue("Fitting model {i} of {n_ns}, node-split: ",
                          as.character(nodesplit$trt2[i]),
                          " vs. ",
                          as.character(nodesplit$trt1[i])))

        ns_arglist$nodesplit <- forcats::fct_c(nodesplit$trt1[i], nodesplit$trt2[i])

        ns_fits[[i + nodesplit_include_consistency]] <- do.call(nma, ns_arglist)
      }

      if (nodesplit_include_consistency) {
        inform(glue::glue("Fitting model {n_ns} of {n_ns}, consistency model"))
        nodesplit <- tibble::add_row(nodesplit, .before = 1)

        ns_arglist$consistency <- "consistency"
        ns_arglist$nodesplit <- NULL

        ns_fits[[1]] <- do.call(nma, ns_arglist)
      }

      nodesplit$model <- ns_fits

      # Return a nma_nodesplit_df object
      class(nodesplit) <- c("nma_nodesplit_df", class(nodesplit))
      return(nodesplit)

    } else if (rlang::is_vector(nodesplit, n = 2)) { # Vector giving single comparison to split

      nodesplit <- as.character(nodesplit)

      if (!all(nodesplit %in% lvls_trt))
        abort(sprintf("The `nodesplit` treatment comparison should match two treatments in the network.\nSuitable values are: %s",
                      ifelse(length(lvls_trt) <= 5,
                             paste0(lvls_trt, collapse = ", "),
                             paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

      if (nodesplit[1] == nodesplit[2])
        abort("`nodesplit` comparison cannot be the same treatment against itself.")

      # Check valid nodesplit - must have both direct and indirect evidence
      if (!has_direct(network, nodesplit[1], nodesplit[2])) {
        abort(glue::glue("Cannot node-split the {nodesplit[1]} vs. {nodesplit[2]} comparison, no direct evidence."))
      }
      if (!has_indirect(network, nodesplit[1], nodesplit[2])) {
        abort(glue::glue("Cannot node-split the {nodesplit[1]} vs. {nodesplit[2]} comparison, no independent indirect evidence."))
      }

      # Store comparison as factor, in increasing order (i.e. trt1 < trt2)
      nodesplit <- factor(nodesplit, levels = lvls_trt)
      if (as.numeric(nodesplit[1]) > as.numeric(nodesplit[2])) {
        nodesplit <- rev(nodesplit)
      }

    } else {
      abort("`nodesplit` should either be a length 2 vector or a 2 column data frame, giving the comparison(s) to node-split.")
    }
  }

  if (!is.null(regression) && !rlang::is_formula(regression, lhs = FALSE)) {
    abort("`regression` should be a one-sided formula.")
  }

  if (is.null(network$classes)) {
    if (!missing(class_interactions)) {
      abort(paste("Setting `class_interactions` requires treatment classes to be specified in the network.",
                  "See set_*() argument `trt_class`.", sep = "\n"))
    } else if (!is.null(regression)) {
      inform(paste("Note: No treatment classes specified in network, any interactions in `regression` formula will be separate (independent) for each treatment.",
                   "Use set_*() argument `trt_class` and nma() argument `class_interactions` to change this.", sep = "\n"))
    }
  }
  class_interactions <- rlang::arg_match(class_interactions)
  if (length(class_interactions) > 1) abort("`class_interactions` must be a single string.")

  likelihood <- check_likelihood(likelihood, network$outcome)
  link <- check_link(link, likelihood)

  # When are priors on auxiliary parameters required?
  has_aux <- (likelihood == "normal" && has_ipd(network)) ||
              likelihood %in% c("ordered", "weibull", "gompertz",
                                "weibull-aft", "lognormal", "loglogistic",
                                "gamma", "gengamma", "mspline", "pexp")

  # Are study intercepts present? Not if only contrast data
  has_intercepts <- has_agd_arm(network) || has_ipd(network)

  # Check priors
  check_prior(prior_intercept)
  check_prior(prior_trt)
  check_prior(prior_het)
  check_prior(prior_reg)
  if (!.is_default(prior_aux)) {
    if (likelihood == "gengamma") check_prior(prior_aux, c("sigma", "k"))
    else check_prior(prior_aux)
  }

  prior_het_type <- rlang::arg_match(prior_het_type)

  # Prior defaults
  prior_defaults <- list()
  if (has_intercepts && .is_default(prior_intercept))
    prior_defaults$prior_intercept <- get_prior_call(prior_intercept)
  if (.is_default(prior_trt))
    prior_defaults$prior_trt <- get_prior_call(prior_trt)
  if (trt_effects == "random" && .is_default(prior_het))
    prior_defaults$prior_het <- get_prior_call(prior_het)
  if (class_effects == "exchangeable" && .is_default(prior_class_mean))
    prior_defaults$prior_class_mean <- get_prior_call(prior_class_mean)
  if (class_effects == "exchangeable" && .is_default(prior_class_sd))
    prior_defaults$prior_class_sd <- get_prior_call(prior_class_sd)
  if (!is.null(regression) && !is_only_offset(regression) && .is_default(prior_reg))
    prior_defaults$prior_reg <- get_prior_call(prior_reg)
  if (has_aux && .is_default(prior_aux)) {
    if (likelihood == "normal" && has_ipd(network)) {
      prior_aux <- .default(half_normal(scale = 5))
    } else if (likelihood == "ordered") {
      prior_aux <- .default(flat())
    } else if (likelihood %in% c("weibull", "gompertz", "weibull-aft",
                                 "lognormal", "loglogistic", "gamma")) {
      prior_aux <- .default(half_normal(scale = 10))
    } else if (likelihood == "gengamma") {
      prior_aux <- .default(list(sigma = half_normal(scale = 10),
                                 k = half_normal(scale = 10)))
    } else if (likelihood %in% c("mspline", "pexp")) {
      prior_aux <- .default(half_normal(scale = 1))
    }
    prior_defaults$prior_aux <- get_prior_call(prior_aux)
  }
  if (has_aux && !is.null(aux_regression) && .is_default(prior_aux_reg)) {
    if (likelihood %in% c("mspline", "pexp")) {
      prior_aux_reg <- .default(half_normal(scale = 1))
    } else if (likelihood %in% valid_lhood$survival) {
      prior_aux_reg <- .default(normal(scale = 10))
    }
    prior_defaults$prior_aux_reg <- get_prior_call(prior_aux_reg)
  }

  # Warn where default priors are used
  if (!rlang::is_empty(prior_defaults)) {
    warn(glue::glue(
      "Prior distributions were left at default values:",
      paste(paste(names(prior_defaults), prior_defaults, sep = " = "), collapse = "\n"),
      .sep = "\n"
    ))
  }

  # Check other args
  if (!rlang::is_bool(QR)) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_bool(center)) abort("`center` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_scalar_integerish(int_thin) ||
      int_thin < 0) abort("`int_thin` should be an integer >= 0.")
  if (!rlang::is_bool(int_check)) abort("`int_check` should be a logical scalar (TRUE or FALSE).")
  if (int_thin > 0) int_check <- FALSE

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!rlang::is_scalar_double(adapt_delta) ||
      adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Check aux_by / aux_regression combination
  aux_by <- rlang::enquo(aux_by)
  if (!rlang::quo_is_null(aux_by) && !is.null(aux_regression)) {
    abort("Cannot specify both `aux_by` and `aux_regression`.")
  }

  if (!is.null(aux_regression)) {
    if (!likelihood %in% valid_lhood$survival) {
      warn("Ignoring `aux_regression`, only supported for survival likelihoods at present.")
    } else if (!has_aux) {
      warn(glue::glue("Ignoring `aux_regression`, no auxiliary parameters in {likelihood} model."))
    }
  }

  if (!rlang::quo_is_null(aux_by)) {
    if (!likelihood %in% valid_lhood$survival) {
      warn("Ignoring `aux_by`, only supported for survival likelihoods at present.")
    } else if (!has_aux) {
      warn(glue::glue("Ignoring `aux_by`, no auxiliary parameters in {likelihood} model."))
    }
  }

  # Set up aux_regression
  if (!is.null(aux_regression) && likelihood %in% valid_lhood$survival &&
      has_aux &&
      (has_ipd(network) || has_agd_arm(network))) {
    if (!rlang::is_formula(aux_regression, lhs = FALSE)) abort("`aux_regression` should be a one-sided formula.")
    if (!is.null(attr(aux_regression, "offset"))) abort("Offset terms not allowed in `aux_regression`.")

    # Remove main effect of study if specified, and make sure treatment handled properly (first, no intercept for M-spline model with symmetric RW prior)
    if (".trt" %in% colnames(attr(terms(aux_regression), "factor"))) {
      aux_regression <- update.formula(aux_regression, ~. -.study -.trt +1)
      if (likelihood %in% c("mspline", "pexp")) aux_regression <- update.formula(aux_regression, ~.trt + . -1)
      else aux_regression <- update.formula(aux_regression, ~.trt + . +1)
    } else {
      aux_regression <- update.formula(aux_regression, ~. -.study +1)
    }
    has_aux_regression <- TRUE
  } else {
    has_aux_regression <- FALSE
    aux_regression <- NULL
  }

  # Set up aux_by
  if (likelihood %in% valid_lhood$survival &&
      has_aux &&
      (has_ipd(network) || has_agd_arm(network))) {

    has_aux_by <- TRUE
    if (rlang::quo_is_null(aux_by)) aux_by <- ".study"

    aux_dat <- dplyr::bind_rows(if (has_ipd(network)) dplyr::select(network$ipd, -".Surv") else NULL,
                                if (has_agd_arm(network)) {
                                  if (inherits(network, "mlnmr_data")) {
                                    .unnest_integration(network$agd_arm) %>% dplyr::select(-".Surv")
                                  } else {
                                    dplyr::select(network$agd_arm, -".Surv")
                                  }
                                } else NULL)

    # Check specs and translate into string column names
    aux_by <- colnames(get_aux_by_data(aux_dat, by = aux_by))
  } else {
    has_aux_by <- FALSE
    aux_by <- NULL
  }


  # Use numerical integration? TRUE if class mlnmr_data and regression is not NULL
  # (Avoids unnecessary use of integration points if regression formula not specified)
  use_int <- inherits(network, "mlnmr_data") && (!is.null(regression) || aux_needs_integration(aux_regression = aux_regression, aux_by = aux_by))

  # Number of numerical integration points
  # Set to 1 if no numerical integration, so that regression on summary data is possible
  n_int <- if (use_int) network$n_int else 1

  # Warn if combining AgD and IPD in a meta-regression without using integration
  if (!is.null(regression) && !inherits(network, "mlnmr_data") && has_ipd(network) &&
      (has_agd_arm(network) || has_agd_contrast(network))) {
    warn(glue::glue("No integration points available, using naive plug-in model at aggregate level.\n",
                    "Use `add_integration()` to add integration points to the network."))
  }

  # Notify if default reference treatment is used
  if (.is_default(network$treatments))
    inform(glue::glue('Note: Setting "{levels(network$treatments)[1]}" as the network reference treatment.'))
  # Notify if network is disconnected
  if (!is_network_connected(network))
    inform("Note: Network is disconnected. See ?is_network_connected for more details.")

  # Get data for design matrices and outcomes
  if (has_ipd(network)) {
    dat_ipd <- network$ipd

    # Only take necessary columns
    dat_ipd <- get_model_data_columns(dat_ipd,
                                      regression = regression,
                                      aux_regression = aux_regression,
                                      label = "IPD",
                                      keep = if (has_aux_by) aux_by else NULL)

    y_ipd <- get_outcome_variables(network$ipd, network$outcome$ipd)
  } else {
    dat_ipd <- tibble::tibble()
    y_ipd <- NULL
  }

  if (has_agd_arm(network)) {
    dat_agd_arm <- network$agd_arm

    y_agd_arm <- get_outcome_variables(dat_agd_arm, network$outcome$agd_arm)

    # Unnest survival data
    if (network$outcome$agd_arm == "survival") {
      y_agd_arm <- tidyr::unnest(y_agd_arm, cols = ".Surv")
    }

    # Set up integration variables if present
    if (use_int) {
      if (network$outcome$agd_arm == "survival") {
        idat_agd_arm <- dat_agd_arm %>%
          # Drop duplicated names in outer dataset from .data_orig before unnesting
          dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(dat_agd_arm))))) %>%
          # Unnest - should have n_int contiguous rows for each survival time
          tidyr::unnest(cols = c(".Surv", ".data_orig")) %>%
          .unnest_integration()
      } else {
        idat_agd_arm <- .unnest_integration(dat_agd_arm)
      }
    } else {
      if (network$outcome$agd_arm == "survival") {
        idat_agd_arm <- dat_agd_arm %>%
          # Drop duplicated names in outer dataset from .data_orig before unnesting
          dplyr::mutate(.data_orig = purrr::map(.data$.data_orig, ~ dplyr::select(., -dplyr::any_of(names(dat_agd_arm))))) %>%
          # Unnest - should have one row for each survival time
          tidyr::unnest(cols = c(".Surv", ".data_orig"))
      } else {
        idat_agd_arm <- dat_agd_arm
      }
    }

    # Only take necessary columns
    idat_agd_arm <- get_model_data_columns(idat_agd_arm,
                                           regression = regression,
                                           aux_regression = aux_regression,
                                           label = "AgD (arm-based)",
                                           keep = if (has_aux_by) aux_by else NULL)

  } else {
    dat_agd_arm <- idat_agd_arm <- tibble::tibble()
    y_agd_arm <- NULL
  }

  if (has_agd_contrast(network)) {
    dat_agd_contrast <- network$agd_contrast

    y_agd_contrast <- get_outcome_variables(dat_agd_contrast, network$outcome$agd_contrast)

    # Set up integration variables if present
    if (use_int) {
      idat_agd_contrast <-  .unnest_integration(dat_agd_contrast)
    } else {
      idat_agd_contrast <- dat_agd_contrast
    }

    # Only take necessary columns
    idat_agd_contrast <- get_model_data_columns(idat_agd_contrast,
                                                regression = regression,
                                                label = "AgD (contrast-based)")

    # Get covariance structure for relative effects, using .se on baseline arm
    Sigma_agd_contrast <- make_Sigma(dat_agd_contrast)

    # Split into baseline and non-baseline arms
    dat_agd_contrast_bl <- dplyr::filter(dat_agd_contrast, is.na(.data$.y))
    idat_agd_contrast_bl <- dplyr::filter(idat_agd_contrast, is.na(.data$.y))
    dat_agd_contrast_nonbl <- dplyr::filter(dat_agd_contrast, !is.na(.data$.y))
    idat_agd_contrast_nonbl <- dplyr::filter(idat_agd_contrast, !is.na(.data$.y))
    y_agd_contrast <- dplyr::filter(y_agd_contrast, !is.na(.data$.y))
  } else {
    dat_agd_contrast <- idat_agd_contrast <-
      dat_agd_contrast_bl <- idat_agd_contrast_bl <-
      dat_agd_contrast_nonbl <- idat_agd_contrast_nonbl <- tibble::tibble()
    y_agd_contrast <- NULL
    Sigma_agd_contrast <- NULL
  }

  # Combine
  idat_all <- dplyr::bind_rows(dat_ipd, idat_agd_arm, idat_agd_contrast_nonbl)
  idat_all_plus_bl <- dplyr::bind_rows(dat_ipd, idat_agd_arm, idat_agd_contrast)

  # Get sample sizes for centering
  if (((!is.null(regression) && !is_only_offset(regression)) || has_aux_regression) && center) {
    # Check that required variables are present in each data set, and non-missing
    if (!is.null(regression)) {
      check_regression_data(regression,
                            dat_ipd = dat_ipd,
                            dat_agd_arm = idat_agd_arm,
                            dat_agd_contrast = idat_agd_contrast)
    }

    if (has_aux_regression) {
      check_regression_data(aux_regression,
                            dat_ipd = dat_ipd,
                            dat_agd_arm = idat_agd_arm,
                            dat_agd_contrast = idat_agd_contrast)
    }

    # If IPD or IPD+AgD use weighted means for centering, otherwise with only AgD use raw mean

    if (has_ipd(network) && (has_agd_arm(network) || has_agd_contrast(network)) && !has_agd_sample_size(network))
      abort(paste("AgD study sample sizes not specified in network, cannot calculate centering values.",
                  "Specify `sample_size` in set_agd_*(), or set center = FALSE.", sep = "\n"))

    if (has_agd_arm(network)) {
      if (has_ipd(network)) {
        N_agd_arm <- network$agd_arm[[".sample_size"]]
      } else {
        N_agd_arm <- rep(n_int, nrow(network$agd_arm))
      }
    } else {
      N_agd_arm <- NULL
    }

    if (has_agd_contrast(network)) {
      if (has_ipd(network)) {
        N_agd_contrast <- network$agd_contrast[[".sample_size"]]
      } else {
        N_agd_contrast <- rep(n_int, nrow(network$agd_contrast))
      }
    } else {
      N_agd_contrast <- NULL
    }

    # Apply weights across integration points
    if (has_agd_arm(network) && network$outcome$agd_arm == "survival") {
      wts <- c(rep(1, nrow(dat_ipd)),
               rep(1 / n_int, nrow(idat_agd_arm)),
               rep(N_agd_contrast / n_int, each = n_int))
    } else {
      wts <- c(rep(1, nrow(dat_ipd)),
               rep(N_agd_arm / n_int, each = n_int),
               rep(N_agd_contrast / n_int, each = n_int))
    }

    # Center numeric columns used in regression model

    if (!is.null(regression)) {
      reg_names <- all.vars(regression)

      # Ignore any variable(s) used as offset(s)
      reg_terms <- terms(regression)

      if (length(attr(reg_terms, "offset"))) {
        off_terms <- rownames(attr(reg_terms, "factors"))[attr(reg_terms, "offset")]
        off_names <- all.vars(as.formula(paste("~", off_terms, sep = "+")))
        reg_names <- setdiff(reg_names, off_names)
      }
    } else {
      reg_names <- character()
    }

    if (has_aux_regression) {
      reg_names <- unique(c(reg_names, all.vars(aux_regression)))
    }

    reg_numeric <- purrr::map_lgl(idat_all[, reg_names], is.numeric)

    # Take weighted mean of all rows (including baseline rows for contrast data)
    if (any(reg_numeric)) {
      xbar <- purrr::map_dbl(idat_all_plus_bl[, reg_names[reg_numeric]], weighted.mean, w = wts)
    } else {
      xbar <- NULL
    }

  } else {
    xbar <- NULL
  }

  # Make NMA formula
  nma_formula <- make_nma_formula(regression,
                                  consistency = consistency,
                                  classes = !is.null(network$classes),
                                  class_interactions = class_interactions)

  # Construct model matrix
  X_list <- make_nma_model_matrix(nma_formula = nma_formula,
                                  dat_ipd = dat_ipd,
                                  dat_agd_arm = idat_agd_arm,
                                  dat_agd_contrast = idat_agd_contrast,
                                  agd_contrast_bl = is.na(idat_agd_contrast$.y),
                                  xbar = xbar,
                                  consistency = consistency,
                                  nodesplit = nodesplit,
                                  classes = !is.null(network$classes))

  X_ipd <- X_list$X_ipd
  X_agd_arm <- X_list$X_agd_arm
  X_agd_contrast <- X_list$X_agd_contrast

  offset_ipd <- X_list$offset_ipd
  offset_agd_arm <- X_list$offset_agd_arm
  offset_agd_contrast <- X_list$offset_agd_contrast

  # Construct RE correlation matrix
  if (trt_effects == "random") {

    # Get study/treatment data
    if (has_ipd(network)) {
      tdat_ipd_arm <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt)
    } else {
      tdat_ipd_arm <- tibble::tibble()
    }

    if (has_agd_arm(network)) {
      tdat_agd_arm <- dplyr::select(dat_agd_arm, ".study", ".trt")
    } else {
      tdat_agd_arm <- tibble::tibble()
    }

    if (has_agd_contrast(network)) {
      tdat_agd_contrast_nonbl <- dplyr::select(dat_agd_contrast_nonbl, ".study", ".trt")
    } else {
      tdat_agd_contrast_nonbl <- tibble::tibble()
    }

    tdat_all <- dplyr::bind_rows(tdat_ipd_arm, tdat_agd_arm, tdat_agd_contrast_nonbl)

    contr <- rep(c(FALSE, FALSE, TRUE),
                 times = c(nrow(tdat_ipd_arm), nrow(tdat_agd_arm), nrow(tdat_agd_contrast_nonbl)))

    if (consistency %in% c("consistency", "nodesplit")) {
      .RE_cor <- RE_cor(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "reftrt")
      .which_RE <- which_RE(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "reftrt")
    } else if (consistency == "ume") {
      .RE_cor <- RE_cor(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "blshift")
      .which_RE <- which_RE(tdat_all$.study, tdat_all$.trt, contrast = contr, type = "blshift")
    } else {
      abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
    }
  } else {
    .RE_cor <- NULL
    .which_RE <- NULL
  }

  # Set up spline basis for mspline and pexp models
  if (likelihood %in% c("mspline", "pexp") && (has_ipd(network) || has_agd_arm(network))) {
    require_pkg("splines2")

    survdat <-
      if (!has_agd_arm(network)) {
        dplyr::tibble(.Surv = y_ipd$.Surv,
                      .study = forcats::fct_drop(dat_ipd$.study),
                      observed = .data$.Surv[, "status"]  == 1)
      } else if (!has_ipd(network)) {
        dplyr::tibble(.Surv = y_agd_arm$.Surv,
                      .study = forcats::fct_drop(rep(dat_agd_arm$.study, times = dat_agd_arm$.sample_size)),
                      observed = .data$.Surv[, "status"]  == 1)
      } else {
        dplyr::tibble(.Surv = c(y_ipd$.Surv, y_agd_arm$.Surv),
                      .study = forcats::fct_drop(forcats::fct_c(dat_ipd$.study, rep(dat_agd_arm$.study, times = dat_agd_arm$.sample_size))),
                      observed = .data$.Surv[, "status"]  == 1)
      }

    survdat <- dplyr::mutate(survdat, !!! get_Surv_data(survdat$.Surv))


    if (!is.null(mspline_basis)) {

      if (!is.list(mspline_basis) || any(!purrr::map_lgl(mspline_basis, is_mspline)))
        abort("`mspline_basis` must be a named list of M-spline bases created using splines2::mSpline().")

      missing_names <- setdiff(levels(survdat$.study), names(mspline_basis))
      if (length(missing_names) > 0)
        abort(glue::glue("`mspline_basis` must be a named list of M-spline basis for each study.\n",
                         "Missing {if (length(missing_names) > 1) 'bases' else 'basis'} for stud{if (length(missing_names) > 1) 'ies' else 'y'} ",
                         glue::glue_collapse(glue::double_quote(missing_names), sep = ", ", last = " and ", width = 30), "."))

      if (!all(purrr::map_int(mspline_basis, ncol) == ncol(mspline_basis[[1]])))
        abort("Each basis in `mspline_basis` must currently have the same number of knots.")

      if (likelihood == "pexp" && any(purrr::map_lgl(mspline_basis, ~attr(., "degree") != 0)))
        abort("`mspline_basis` for piecewise exponential model must all have degree = 0.")

      basis <- mspline_basis

    } else {

      if (likelihood == "pexp") mspline_degree <- 0

      if (!rlang::is_scalar_integerish(mspline_degree, finite = TRUE) || mspline_degree < 0)
        abort("`mspline_degree` must be a single non-negative integer.")

      if (is.null(knots)) {

        knots <- make_knots(network, n_knots = n_knots, type = if (!is.null(aux_regression)) "quantile_common" else "quantile")

      } else {  # User-provided knots
        # Check required format
        if (!has_aux_regression) {
          if ((!is.list(knots) || any(!purrr::map_lgl(knots, is.numeric))) && !is.vector(knots, mode = "numeric"))
            abort("`knots` must be a named list of numeric vectors giving the knot locations for each study, or a single numeric vector of knot locations to use for all studies.")

          if (is.vector(knots, mode = "numeric")) {
            knots <- rep_len(list(knots), dplyr::n_distinct(survdat$.study))
            names(knots) <- unique(survdat$.study)
          }

          missing_names <- setdiff(levels(survdat$.study), names(knots))
          if (length(missing_names) > 0)
            abort(glue::glue("`knots` must be a named list of numeric vectors giving the knot locations for each study.\n",
                             "Missing knot location vector{if (length(missing_names) > 1) 's' else ''} for stud{if (length(missing_names) > 1) 'ies' else 'y'} ",
                             glue::glue_collapse(glue::double_quote(missing_names), sep = ", ", last = " and ", width = 30), "."))

          if (!all(purrr::map_int(knots, length) == length(knots[[1]])))
            abort("Each vector in `knots` must currently be the same length (each study must have the same number of knots).")

          if (length(knots[[1]]) < 3)
            abort("Each vector in `knots` must be at least length 3 (boundary knots and one internal knot).")
        } else {
          # With aux_regression, only a single vector of knot locations is allowed
          if (!is.vector(knots, mode = "numeric"))
            abort("`knots` must be a single numeric vector of knot locations, shared for all studies, when `aux_regression` is specified.")

          if (length(knots) < 3)
            abort("`knots` must be at least length 3 (boundary knots and one internal knot).")

          knots <- rep_len(list(knots), dplyr::n_distinct(survdat$.study))
          names(knots) <- unique(survdat$.study)
        }
      }

      # Set up basis
      # Only evaluate at first boundary knot for now to save time/memory
      knots <- purrr::map(knots, sort)
      knots <- tibble::as_tibble(knots)
      b_knots <- knots[c(1, nrow(knots)), ]
      i_knots <- knots[-c(1, nrow(knots)), ]
      basis <- purrr::imap(b_knots,
                           ~withCallingHandlers(splines2::mSpline(.x[1],
                                                       knots = i_knots[[.y]],
                                                       Boundary.knots = .x,
                                                       degree = mspline_degree,
                                                       intercept = TRUE),
                                     error = function(e) abort(glue::glue("Could not create spline basis for study {glue::double_quote(.y)}."),
                                                               parent = e),
                                     warning = function(w) {
                                       warn(glue::glue("Warning while creating spline basis for study {glue::double_quote(.y)}."), parent = w)
                                       rlang::cnd_muffle(w)
                                     }
                                    )
                           )

    }

    # Ensure list is in factor order
    basis <- basis[levels(survdat$.study)]

  } else {
    basis <- NULL
  }

  # Set up aux_by design vector
  if (has_aux_by) {
    # Determine whether we need to integrate over the aux regression too
    aux_int <- aux_needs_integration(aux_regression = aux_regression, aux_by = aux_by)

    if (aux_int) {
      aux_dat <- dplyr::bind_rows(dat_ipd, idat_agd_arm)
    } else {
      aux_dat <- dplyr::bind_rows(dat_ipd,
                                  if (has_agd_arm(network)) dplyr::select(tidyr::unnest(dat_agd_arm, cols = ".Surv"), -".Surv") else NULL)
    }
    aux_id <- get_aux_id(aux_dat, aux_by)
  } else {
    aux_id <- integer()
    aux_int <- FALSE
  }

  # Set up aux_regression design matrix
  if (has_aux_regression) {
    X_aux <- make_nma_model_matrix(aux_regression,
                                   dat_ipd = aux_dat,
                                   xbar = xbar)$X_ipd

    # Group common rows of X_aux for efficiency if possible
    if (aux_int && use_int && has_agd_arm(network)) {
      aux_group <- 1:nrow(X_aux)
    } else {
      X_aux_dat <- as.data.frame(X_aux)
      if (!rlang::has_name(X_aux_dat, ".study")) X_aux_dat <- dplyr::mutate(X_aux_dat, aux_dat$.study, .before = 0)  # Always group within study
      aux_group <- dplyr::group_indices(dplyr::group_by_all(X_aux_dat))
    }
  } else {
    X_aux <- NULL
    aux_group <- aux_id
  }

if (class_effects == "exchangeable") {
  # Create class design vector for class means
  class_mean_design <- which_CE(network$classes, class_sd)

  # Create class design vector for class SDs
  if (is.list(class_sd)) {
    class_sd_design <- which_CE(forcats::fct_collapse(network$classes, !!!class_sd), class_sd)
  } else if (class_sd == "common") {
    class_sd_design <- list(
      # Change non-zero class IDs to 1
      id = pmin(class_mean_design$id, 1),
      # Set common class label
      label = "All Classes"
    )
  } else if (class_sd == "independent") {
    class_sd_design <- class_mean_design
  }
}
  # Fit using nma.fit
  stanfit <- nma.fit(ipd_x = X_ipd, ipd_y = y_ipd,
    agd_arm_x = X_agd_arm, agd_arm_y = y_agd_arm,
    agd_contrast_x = X_agd_contrast, agd_contrast_y = y_agd_contrast,
    agd_contrast_Sigma = Sigma_agd_contrast,
    n_int = n_int,
    ipd_offset = offset_ipd,
    agd_arm_offset = offset_agd_arm,
    agd_contrast_offset = offset_agd_contrast,
    trt_effects = trt_effects,
    RE_cor = .RE_cor,
    which_RE = .which_RE,
    class_effects = class_effects,
    which_CE = if (class_effects == "exchangeable") class_mean_design$id else NULL,
    which_CE_sd = if (class_effects == "exchangeable") class_sd_design$id else NULL,
    likelihood = likelihood,
    link = link,
    consistency = consistency,
    ...,
    prior_intercept = prior_intercept,
    prior_trt = prior_trt,
    prior_het = prior_het,
    prior_het_type = prior_het_type,
    prior_reg = prior_reg,
    prior_aux = prior_aux,
    prior_aux_reg = prior_aux_reg,
    prior_class_mean = prior_class_mean,
    prior_class_sd = prior_class_sd,
    aux_id = aux_id,
    aux_group = aux_group,
    X_aux = X_aux,
    QR = QR,
    adapt_delta = adapt_delta,
    int_thin = int_thin,
    int_check = int_check,
    basis = basis)

  # Make readable parameter names for generated quantities
  fnames_oi <- stanfit@sim$fnames_oi

  # Labels for fitted values
  ipd_data_labels <- if (has_ipd(network)) make_data_labels(dat_ipd$.study, dat_ipd$.trt) else NULL

  if (has_agd_arm(network)) {
    if (network$outcome$agd_arm != "survival") {
      agd_arm_data_labels <-  make_data_labels(dat_agd_arm$.study, dat_agd_arm$.trt)
    } else {
      agd_arm_data_labels <-  make_data_labels(rep(dat_agd_arm$.study, times = dat_agd_arm$.sample_size),
                                               rep(dat_agd_arm$.trt, times = dat_agd_arm$.sample_size))
    }
  } else {
    agd_arm_data_labels <- NULL
  }

  if (has_agd_contrast(network)) {
    dat_agd_contrast_nonbl <-
      dplyr::left_join(dat_agd_contrast_nonbl,
                       dplyr::transmute(dat_agd_contrast_bl, .data$.study, .trt_b = .data$.trt),
                       by = ".study")

    agd_contrast_data_labels <- make_data_labels(dat_agd_contrast_nonbl$.study,
                                                 dat_agd_contrast_nonbl$.trt,
                                                 dat_agd_contrast_nonbl$.trt_b)
  } else {
    agd_contrast_data_labels <- NULL
  }

  if (has_ipd(network))
    fnames_oi[grepl("^fitted_ipd\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_ipd[", ipd_data_labels, "]")

  if (has_agd_arm(network))
    fnames_oi[grepl("^fitted_agd_arm\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_agd_arm[", agd_arm_data_labels, "]")

  if (has_agd_contrast(network))
    fnames_oi[grepl("^fitted_agd_contrast\\[[0-9]+\\]$", fnames_oi)] <- paste0("fitted_agd_contrast[", agd_contrast_data_labels, "]")

  # Labels for RE deltas
  if (trt_effects == "random") {

    if (consistency == "ume") { # Baseline shift models (currently only UME)
      if (has_ipd(network)) {
        ipd_arms <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt) %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::mutate(.trt_b = sort(.data$.trt)[1])
        ipd_delta_labels <- make_data_labels(ipd_arms$.study, ipd_arms$.trt, ipd_arms$.trt_b)
      } else {
        ipd_delta_labels <- NULL
      }

      if (has_agd_arm(network)) {
        agd_arm_trt_b <- dat_agd_arm %>%
          dplyr::group_by(.data$.study) %>%
          dplyr::mutate(.trt_b = sort(.data$.trt)[1]) %>%
          dplyr::pull(.data$.trt_b)
        agd_arm_delta_labels <- make_data_labels(dat_agd_arm$.study, dat_agd_arm$.trt, agd_arm_trt_b)
      } else {
        agd_arm_delta_labels <- NULL
      }

    } else { # Reference treatment models
      if (has_ipd(network)) {
        ipd_arms <- dplyr::distinct(dat_ipd, .data$.study, .data$.trt)
        ipd_delta_labels <- make_data_labels(ipd_arms$.study, ipd_arms$.trt)
      } else {
        ipd_delta_labels <- NULL
      }
      agd_arm_delta_labels <- if (has_agd_arm(network)) make_data_labels(dat_agd_arm$.study, dat_agd_arm$.trt) else NULL
    }

    agd_contrast_delta_labels <- agd_contrast_data_labels

    delta_labels <- c(ipd_delta_labels,
                      agd_arm_delta_labels,
                      agd_contrast_delta_labels)[.which_RE > 0]

    fnames_oi[grepl("^delta\\[[0-9]+\\]$", fnames_oi)] <- paste0("delta[", delta_labels, "]")
  }

  # Labels for log_lik, resdev (only one entry per AgD contrast study)
  dev_labels <- c(ipd_data_labels,
                  agd_arm_data_labels,
                  if (has_agd_contrast(network)) as.character(dat_agd_contrast_bl$.study) else NULL)

  fnames_oi[grepl("^log_lik\\[[0-9]+\\]$", fnames_oi)] <- paste0("log_lik[", dev_labels, "]")
  fnames_oi[grepl("^resdev\\[[0-9]+\\]$", fnames_oi)] <- paste0("resdev[", dev_labels, "]")

  # Labels for cumulative integration points
  if (inherits(network, "mlnmr_data") && (has_agd_arm(network) || has_agd_contrast(network)) && int_thin > 0) {
    n_int_thin <- n_int %/% int_thin

    if (has_agd_arm(network)) {
      agd_arm_cumint_labels <- rep(agd_arm_data_labels, each = n_int_thin)
      agd_arm_cumint_labels <- paste0(agd_arm_cumint_labels, ", ", rep_len(1:n_int_thin * int_thin, length.out = length(agd_arm_cumint_labels)))

      if (likelihood == "ordered") {
        l_cat <- colnames(y_agd_arm$.r)
        agd_arm_cumint_labels <- paste0(rep(agd_arm_cumint_labels, times = length(l_cat)), ", ", rep(l_cat, each = length(agd_arm_cumint_labels)))
      }

      fnames_oi[grepl("^theta_bar_cum_agd_arm\\[", fnames_oi)] <- paste0("theta_bar_cum_agd_arm[", agd_arm_cumint_labels, "]")
      if (likelihood %in% c("bernoulli2", "binomial2"))
        fnames_oi[grepl("^theta2_bar_cum\\[[0-9]+\\]$", fnames_oi)] <- paste0("theta2_bar_cum[", agd_arm_cumint_labels, "]")
    }

    if (has_agd_contrast(network)) {
      agd_contrast_cumint_labels <- rep(agd_contrast_data_labels, each = n_int_thin)
      agd_contrast_cumint_labels <- paste0(agd_contrast_cumint_labels, ", ", rep_len(1:n_int_thin * int_thin, length.out = length(agd_contrast_cumint_labels)))

      fnames_oi[grepl("^theta_bar_cum_agd_contrast\\[[0-9]+\\]$", fnames_oi)] <- paste0("theta_bar_cum_agd_contrast[", agd_contrast_cumint_labels, "]")
    }

  }

  # Labels for survival aux pars
  if (likelihood %in% valid_lhood$survival && has_aux) {
    aux_labels <- get_aux_labels(aux_dat, by = aux_by)

    fnames_oi[grepl("^shape\\[[0-9]+\\]$", fnames_oi)] <- paste0("shape[", aux_labels, "]")
    fnames_oi[grepl("^sdlog\\[[0-9]+\\]$", fnames_oi)] <- paste0("sdlog[", aux_labels, "]")
    fnames_oi[grepl("^sigma\\[[0-9]+\\]$", fnames_oi)] <- paste0("sigma[", aux_labels, "]")
    fnames_oi[grepl("^k\\[[0-9]+\\]$", fnames_oi)] <- paste0("k[", aux_labels, "]")

    if (likelihood %in% c("mspline", "pexp")) {
      # Number of spline coefficients
      n_scoef <- ncol(basis[[1]])

      fnames_oi[grepl("^scoef\\[[0-9]+,[0-9]+\\]$", fnames_oi)] <-
        paste0("scoef[", rep(aux_labels, times = n_scoef), ", ", rep(1:n_scoef, each = length(aux_labels)), "]")
    }
  }

  if (class_effects == "exchangeable"){
    # Label class_mean parameters
    fnames_oi[grepl("^class_mean\\[[0-9]+\\]$", fnames_oi)] <- paste0("class_mean[", class_mean_design$label, "]")

    # Label class_sd parameters
    fnames_oi[grepl("^class_sd\\[[0-9]+\\]$", fnames_oi)] <- paste0("class_sd[", class_sd_design$label, "]")
    network$class_sd <- class_sd_design$label
}
  stanfit@sim$fnames_oi <- fnames_oi

  # Create stan_nma object
  out <- list(network = network,
              stanfit = stanfit,
              trt_effects = trt_effects,
              consistency = consistency,
              regression = regression,
              aux_regression = aux_regression,
              class_interactions = if (!is.null(regression) && !is.null(network$classes)) class_interactions else NULL,
              class_effects = class_effects,
              class_sd = if (class_effects == "exchangeable") class_sd else NULL,
              xbar = xbar,
              likelihood = likelihood,
              link = link,
              aux_by = if (has_aux_by) colnames(get_aux_by_data(aux_dat, by = aux_by)) else NULL,
              priors = list(prior_intercept = if (has_intercepts) prior_intercept else NULL,
                            prior_trt = prior_trt,
                            prior_class_mean = if (class_effects == "exchangeable") prior_class_mean else NULL,
                            prior_class_sd = if (class_effects == "exchangeable") prior_class_sd else NULL,
                            prior_het = if (trt_effects == "random") prior_het else NULL,
                            prior_het_type = if (trt_effects == "random") prior_het_type else NULL,
                            prior_reg = if (!is.null(regression) && !is_only_offset(regression)) prior_reg else NULL,
                            prior_aux = if (has_aux) prior_aux else NULL,
                            prior_aux_reg = if (has_aux_regression) prior_aux_reg else NULL))

  if (likelihood %in% c("mspline", "pexp")) out$basis <- basis

  if (inherits(network, "mlnmr_data")) class(out) <- c("stan_mlnmr", "stan_nma")
  else class(out) <- "stan_nma"

  if (likelihood %in% valid_lhood$survival) class(out) <- c("stan_nma_surv", class(out))

  if (consistency == "nodesplit" && !is.data.frame(nodesplit)) {
    class(out) <- c("nma_nodesplit", class(out))
    out$nodesplit <- nodesplit
  }

  return(out)
}


#' @param ipd_x Design matrix for IPD studies
#' @param ipd_y Outcome data frame for IPD studies
#' @param agd_arm_x  Design matrix for AgD studies (arm-based)
#' @param agd_arm_y  Outcome data frame for AgD studies (arm-based)
#' @param agd_contrast_x  Design matrix for AgD studies (contrast-based)
#' @param agd_contrast_y  Outcome data frame for AgD studies (contrast-based)
#' @param agd_contrast_Sigma List of covariance matrices for contrast-based data
#' @param n_int Number of numerical integration points used
#' @param ipd_offset Vector of offset values for IPD
#' @param agd_arm_offset Vector of offset values for AgD (arm-based)
#' @param agd_contrast_offset Vector of offset values for AgD (contrast-based)
#' @param RE_cor Random effects correlation matrix, when `trt_effects = "random"`
#' @param which_RE Random effects design vector, when `trt_effects = "random"`
#' @param which_CE Class effects means design vector (0 = no class)
#' @param which_CE_sd Class effects SDs design vector (0 = no class)
#' @param basis Spline basis for `mspline` and `pexp` models
#'
#' @noRd
nma.fit <- function(ipd_x, ipd_y,
                    agd_arm_x, agd_arm_y,
                    agd_contrast_x, agd_contrast_y, agd_contrast_Sigma,
                    n_int,
                    ipd_offset = NULL, agd_arm_offset = NULL, agd_contrast_offset = NULL,
                    trt_effects = c("fixed", "random"),
                    RE_cor = NULL,
                    which_RE = NULL,
                    class_effects = c("independent", "exchangeable", "common"),
                    which_CE = NULL,
                    which_CE_sd = NULL,
                    likelihood = NULL,
                    link = NULL,
                    consistency = c("consistency", "ume", "nodesplit"),
                    ...,
                    prior_intercept,
                    prior_trt,
                    prior_het,
                    prior_het_type = c("sd", "var", "prec"),
                    prior_reg,
                    prior_aux,
                    prior_aux_reg,
                    prior_class_mean,
                    prior_class_sd,
                    aux_id = integer(),
                    aux_group = integer(),
                    X_aux = NULL,
                    QR = FALSE,
                    adapt_delta = NULL,
                    int_thin = 0,
                    int_check = TRUE,
                    basis) {

  if (missing(ipd_x)) ipd_x <- NULL
  if (missing(ipd_y)) ipd_y <- NULL
  if (missing(agd_arm_x)) agd_arm_x <- NULL
  if (missing(agd_arm_y)) agd_arm_y <- NULL
  if (missing(agd_contrast_x)) agd_contrast_x <- NULL
  if (missing(agd_contrast_y)) agd_contrast_y <- NULL
  if (missing(agd_contrast_Sigma)) agd_contrast_Sigma <- NULL

  # Check available x and y
  if (xor(is.null(ipd_x), is.null(ipd_y)))
    abort("`ipd_x` and `ipd_y` should both be present or both NULL.")
  if (xor(is.null(agd_arm_x), is.null(agd_arm_y)))
    abort("`agd_arm_x` and `agd_arm_y` should both be present or both NULL.")
  if (xor(is.null(agd_contrast_x), is.null(agd_contrast_y)) || xor(is.null(agd_contrast_x), is.null(agd_contrast_Sigma)))
    abort("`agd_contrast_x`, `agd_contrast_y`, `agd_contrast_Sigma` should all be present or all NULL.")

  has_ipd <- !is.null(ipd_x) && !is.null(ipd_y)
  has_agd_arm <- !is.null(agd_arm_x) && !is.null(agd_arm_y)
  has_agd_contrast <-  !is.null(agd_contrast_x) && !is.null(agd_contrast_y) && !is.null(agd_contrast_Sigma)

  # Ignore n_int if no AgD
  if (!has_agd_arm && !has_agd_contrast) n_int <- 1
  # Check n_int
  if (!rlang::is_scalar_integerish(n_int) ||
      n_int < 1) abort("`n_int` should be an integer >= 1.")

  # Check design matrices, outcomes
  if (has_ipd) {
    if (!is.matrix(ipd_x) || !is.numeric(ipd_x))
      abort("`ipd_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(ipd_y, is.numeric)))
      abort("`ipd_y` should be numeric outcome data.")
    if (nrow(ipd_x) != nrow(ipd_y))
      abort("Number of rows in `ipd_x` and `ipd_y` do not match.")
  }
  if (has_agd_arm) {
    if (!is.matrix(agd_arm_x) || !is.numeric(agd_arm_x))
      abort("`agd_arm_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_arm_y, is.numeric)))
      abort("`agd_arm_y` should be numeric outcome data.")
    if (nrow(agd_arm_x) != nrow(agd_arm_y) * n_int)
      abort("Number of rows in `agd_arm_x`, `agd_arm_y`, and `n_int` do not match.")
  }
  if (has_agd_contrast) {
    if (!is.matrix(agd_contrast_x) || !is.numeric(agd_contrast_x))
      abort("`agd_contrast_x` should be a numeric matrix.")
    if (any(!purrr::map_lgl(agd_contrast_y, is.numeric)))
      abort("`agd_contrast_y` should be numeric outcome data.")
    if (nrow(agd_contrast_x) != nrow(agd_contrast_y) * n_int)
      abort("Number of rows in `agd_contrast_x`, `agd_contrast_y`, and `n_int` do not match.")
    if (!is.list(agd_contrast_Sigma) || any(purrr::map_lgl(agd_contrast_Sigma, ~!is.numeric(.))))
      abort("`agd_contrast_Sigma` should be a list of covariance matrices, of length equal to the number of AgD (contrast-based) studies.")
  }

  # Check offsets
  has_ipd_offset <- !is.null(ipd_offset)
  has_agd_arm_offset <- !is.null(agd_arm_offset)
  has_agd_contrast_offset <- !is.null(agd_contrast_offset)

  if (has_ipd_offset && !has_ipd)
    abort("`ipd_offset` given but not `ipd_x` and `ipd_y`.")
  if (has_agd_arm_offset && !has_agd_arm)
    abort("`agd_arm_offset` given but not `agd_arm_x` and `agd_arm_y`.")
  if (has_agd_contrast_offset && !has_agd_contrast)
    abort("`agd_contrast_offset` given but not `agd_contrast_x` and `agd_contrast_y`.")

  has_offsets <- any(has_ipd_offset, has_agd_arm_offset, has_agd_contrast_offset)
  if (has_offsets &&
      !all(has_ipd_offset[has_ipd],
           has_agd_arm_offset[has_agd_arm],
           has_agd_contrast_offset[has_agd_contrast]))
    abort("Offsets provided for some data sources but not all.")

  if (has_ipd_offset && (!is.numeric(ipd_offset) || length(ipd_offset) != nrow(ipd_x)))
    abort("`ipd_offset` should be a numeric vector with length matching `ipd_x`, `ipd_y`")
  if (has_agd_arm_offset && (!is.numeric(agd_arm_offset) || length(agd_arm_offset) != nrow(agd_arm_x)))
    abort("`agd_arm_offset` should be a numeric vector with length matching `agd_arm_x`, `agd_arm_y`")
  if (has_agd_contrast_offset && (!is.numeric(agd_contrast_offset) || length(agd_contrast_offset) != nrow(agd_contrast_x)))
    abort("`agd_contrast_offset` should be a numeric vector with length matching `agd_contrast_x`, `agd_contrast_y`")

  # Check matching X column names and dimensions if more than one present
  if ((has_ipd &&
       ((has_agd_arm && !identical(colnames(ipd_x), colnames(agd_arm_x))) ||
        (has_agd_contrast && !identical(colnames(ipd_x), colnames(agd_contrast_x))))) ||
      (has_agd_arm &&
       (has_agd_contrast && !identical(colnames(agd_arm_x), colnames(agd_contrast_x)))))
    abort("Non-matching columns in *_x matrices.")

  # Check model arguments
  trt_effects <- rlang::arg_match(trt_effects)
  if (length(trt_effects) > 1) abort("`trt_effects` must be a single string.")

  # Check class effect arguments
  class_effects <- rlang::arg_match(class_effects)
  if (length(class_effects) > 1) abort("`class_effects` must be a single string.")
if (class_effects == "exchangeable") {
  if (is.null(which_CE) || !rlang::is_integerish(which_CE) || any(which_CE < 0))
    abort("`which_CE` must be an integer design vector for class effects.")
  if (is.null(which_CE_sd) || !rlang::is_integerish(which_CE_sd) || any(which_CE_sd < 0))
    abort("`which_CE_sd` must be an integer design vector for class effect SDs.")
}

  likelihood <- check_likelihood(likelihood)
  link <- check_link(link, likelihood)

  # When are priors on auxiliary parameters required?
  has_aux <- (likelihood == "normal" && has_ipd) ||
    (likelihood %in% c("ordered", setdiff(valid_lhood$survival, c("exponential", "exponential-aft"))) &&
       (has_ipd || has_agd_arm))

  # Check priors
  check_prior(prior_intercept)
  check_prior(prior_trt)
  if (trt_effects == "random") check_prior(prior_het)
  check_prior(prior_reg)
  if (has_aux) {
    if (likelihood == "gengamma") check_prior(prior_aux, c("sigma", "k"))
    else check_prior(prior_aux)

    if (!is.null(X_aux)) check_prior(prior_aux_reg)
  }
  if (class_effects == "exchangeable"){
  check_prior(prior_class_mean)
  check_prior(prior_class_sd)
} else {
  # Dummy class effects priors for non-CE models, not used but requested by Stan data
  prior_class_mean <- normal(0, 1)
  prior_class_sd <- half_normal(1)
}
  prior_het_type <- rlang::arg_match(prior_het_type)

  # Dummy RE prior for FE model, not used but requested by Stan data
  if (trt_effects == "fixed") prior_het <- half_normal(1)

  # Check other args
  if (!rlang::is_bool(QR)) abort("`QR` should be a logical scalar (TRUE or FALSE).")
  if (!rlang::is_scalar_integerish(int_thin) ||
      int_thin < 0) abort("`int_thin` should be an integer >= 0.")
  if (!rlang::is_bool(int_check)) abort("`int_check` should be a logical scalar (TRUE or FALSE).")
  if (int_thin > 0) int_check <- FALSE

  # Set adapt_delta
  if (is.null(adapt_delta)) {
    adapt_delta <- switch(trt_effects, fixed = 0.8, random = 0.95)
  } else if (!is.numeric(adapt_delta) ||
             length(adapt_delta) > 1 ||
             adapt_delta <= 0 || adapt_delta >= 1) abort("`adapt_delta` should be a  numeric value in (0, 1).")

  # Is this a survival outcome?
  is_survival <- likelihood %in% valid_lhood$survival

  # Pull study and treatment details from *_x
  if (has_ipd) x_names <- colnames(ipd_x)
  else if (has_agd_arm) x_names <- colnames(agd_arm_x)
  else if (has_agd_contrast) x_names <- colnames(agd_contrast_x)

  col_study <- grepl("^\\.study[^:]+$", x_names)
  col_trt <- grepl("^(\\.trt|\\.contr)[^:]+$", x_names)
  col_omega <- x_names == ".omegaTRUE"
  col_reg <- !col_study & !col_trt & !col_omega

  n_trt <- sum(col_trt) + 1

  get_study <- function(x) which(x == 1)
  get_trt <- function(x, v = 1) if (any(x == v)) which(x == v) + 1 else 1

  if (has_ipd) {
    ipd_s_t_all <- dplyr::tibble(.study = unname(apply(ipd_x[, col_study, drop = FALSE], 1, get_study)),
                                 .trt = unname(apply(ipd_x[, col_trt, drop = FALSE], 1, get_trt)))
    ipd_s_t <- dplyr::distinct(ipd_s_t_all) %>% dplyr::mutate(.arm = 1:dplyr::n())
    ipd_arm <-  dplyr::left_join(ipd_s_t_all, ipd_s_t, by = c(".study", ".trt")) %>% dplyr::pull(.data$.arm)
    ipd_study <- ipd_s_t$.study
    ipd_trt <- ipd_s_t$.trt
    narm_ipd <- max(ipd_arm)
    ni_ipd <- nrow(ipd_x)
  } else {
    ipd_s_t_all <- dplyr::tibble(.study = integer(), .trt = integer())
    ipd_study <- ipd_trt <- ipd_arm <- numeric()
    ni_ipd <- 0
    narm_ipd <- 0
  }

  if (has_agd_arm) {
    ni_agd_arm <- nrow(agd_arm_y)
    aa1 <- 0:(ni_agd_arm - 1)*n_int + 1
    agd_arm_study <- apply(agd_arm_x[aa1, col_study, drop = FALSE], 1, get_study)
    agd_arm_trt <- apply(agd_arm_x[aa1, col_trt, drop = FALSE], 1, get_trt)

    if (!is_survival) {
      narm_agd_arm <- ni_agd_arm
    } else {
      agd_arm_s_t_all <- dplyr::tibble(.study = agd_arm_study, .trt = agd_arm_trt)
      agd_arm_s_t <- dplyr::distinct(agd_arm_s_t_all) %>% dplyr::mutate(.arm = 1:dplyr::n())
      agd_arm_arm <-  dplyr::left_join(agd_arm_s_t_all, agd_arm_s_t, by = c(".study", ".trt")) %>% dplyr::pull(.data$.arm)
      agd_arm_study <- agd_arm_s_t$.study
      agd_arm_trt <- agd_arm_s_t$.trt
      narm_agd_arm <- max(agd_arm_arm)
    }

  } else {
    agd_arm_s_t_all <- dplyr::tibble(.study = integer(), .trt = integer())
    agd_arm_study <- agd_arm_trt <- agd_arm_arm <- numeric()
    ni_agd_arm <- narm_agd_arm <- 0
  }

  if (has_agd_contrast) {
    ni_agd_contrast <- nrow(agd_contrast_y)
    ac1 <- 0:(ni_agd_contrast - 1)*n_int + 1
    agd_contrast_trt <- apply(agd_contrast_x[ac1, col_trt, drop = FALSE], 1, get_trt)
    agd_contrast_trt_b <- apply(agd_contrast_x[ac1, col_trt, drop = FALSE], 1, get_trt, v = -1)

    # Get number of contrast-based studies from length of Sigma list
    ns_agd_contrast <- length(agd_contrast_Sigma)

    # Construct block-diagonal contrast covariance matrix
    Sigma <- as.matrix(Matrix::bdiag(agd_contrast_Sigma))

    if (nrow(Sigma) != ni_agd_contrast)
      abort("Dimensions of `agd_contrast_Sigma` covariance matrices do not match the contrast-based data.")
  } else {
    agd_contrast_trt <- agd_contrast_trt_b <- numeric()
    Sigma <- matrix(1, 1, 1)
    ni_agd_contrast <- 0
    ns_agd_contrast <- 0
  }

  # Set up random effects
  if (trt_effects == "random") {
    narm <- narm_ipd + narm_agd_arm + ni_agd_contrast
    if (!is.null(which_RE)) {
      if (!rlang::is_integerish(which_RE) ||
          any(which_RE < 0) ||
          is.matrix(which_RE))
        abort("`which_RE` should be an integer vector.")
      if (length(which_RE) != narm)
        abort(glue::glue("Length of `which_RE` does not match the number of arms/contrasts.\n",
                         "Expecting length {narm}, instead length {length(which_RE)}."))
    } else {
      abort("Specify `which_RE` when trt_effects = 'random'.")
    }

    nRE <- sum(which_RE != 0)

    if (!is.null(RE_cor)) {
      if (!is.matrix(RE_cor) || !is.numeric(RE_cor))
        abort("`RE_cor` should be a numeric matrix.")
      if (any(dim(RE_cor) != c(nRE, nRE)))
        abort(glue::glue("Dimensions of `RE_cor` do not match the number of random effects.\n",
                         "Expecting [{nRE} x {nRE}], instead [{nrow(RE_cor)} x {ncol(RE_cor)}]."))
    } else {
      abort("Specify `RE_cor` when trt_effects = 'random'.")
    }
  } else {
    RE_cor <- matrix(1, 1, 1)
    which_RE <- integer()
  }

  # Make full design matrix
  X_all <- rbind(ipd_x, agd_arm_x, agd_contrast_x)

  # Make sure columns of X_all are in correct order (study, trt, omega (if nodesplit), regression terms)
  X_all <- cbind(X_all[, col_study, drop = FALSE],
                 X_all[, col_trt, drop = FALSE],
                 X_all[, col_omega, drop = FALSE],
                 X_all[, col_reg, drop = FALSE])

  # Take thin QR decomposition if QR = TRUE
  if (QR) {
    X_all_qr <- qr(X_all)
    X_all_Q <- qr.Q(X_all_qr) * sqrt(nrow(X_all) - 1)
    X_all_R <- qr.R(X_all_qr)[, sort.list(X_all_qr$pivot)] / sqrt(nrow(X_all) - 1)
    X_all_R_inv <- solve(X_all_R)
  }

  # Handle integration points
  dots <- list(...)
  nchains <- dots$chains
  if (is.null(nchains)) nchains <- 4
  if (nchains < 2) {
    warn("At least 2 chains are required to check integration convergence with int_check = TRUE.")
    int_check <- FALSE
  }
  nint_max <- n_int
  if (int_check && n_int > 1) {
    nint_vec <- c(rep(n_int, ceiling(nchains / 2)), rep(n_int %/% 2, floor(nchains / 2)))
  } else {
    nint_vec <- rep(n_int, nchains)
  }

  # Set common Stan data
  standat <- list(
    # Constants
    ns_ipd = length(unique(ipd_study)),
    ni_ipd = ni_ipd,
    ns_agd_arm = length(unique(agd_arm_study)),
    ni_agd_arm = ni_agd_arm,
    ns_agd_contrast = ns_agd_contrast,
    ni_agd_contrast = ni_agd_contrast,
    nt = n_trt,
    nchains = nchains,
    nint_max = nint_max,
    nint_vec = nint_vec,
    nX = ncol(X_all),
    int_thin = int_thin,
    # Study and treatment details
    narm_ipd = narm_ipd,
    ipd_arm = ipd_arm,
    ipd_trt = ipd_trt,
    narm_agd_arm = narm_agd_arm,
    agd_arm_trt = agd_arm_trt,
    agd_contrast_trt = as.array(agd_contrast_trt),
    agd_contrast_trt_b = as.array(agd_contrast_trt_b),
    agd_contrast_y = if (has_agd_contrast) as.array(agd_contrast_y$.y) else numeric(),
    agd_contrast_Sigma = Sigma,
    # ipd_study = ipd_study,
    # agd_arm_study = agd_arm_study,
    # agd_contrast_study = agd_contrast_study,
    # Random effects
    RE = switch(trt_effects, fixed = 0, random = 1),
    RE_cor = RE_cor,
    which_RE = which_RE,
    # Node splitting
    nodesplit = consistency == "nodesplit",
    # Design matrix or QR decomposition
    QR = QR,
    X = if (QR) X_all_Q else X_all,
    R_inv = if (QR) X_all_R_inv else matrix(0, 0, 0),
    # Offsets
    has_offset = has_offsets,
    offsets = if (has_offsets) as.array(c(ipd_offset, agd_arm_offset, agd_contrast_offset)) else numeric(),
    # Class effects
    which_CE = if (class_effects == "exchangeable") which_CE else numeric(0),
    which_CE_sd = if (class_effects == "exchangeable") which_CE_sd else numeric(0),
    class_effects = ifelse(class_effects == "exchangeable", 1, 0)
    )

  # Add priors
  standat <- purrr::list_modify(standat,
    !!! prior_standat(prior_intercept, "prior_intercept",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_trt, "prior_trt",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_reg, "prior_reg",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_het, "prior_het",
                      valid = c("Normal", "half-Normal", "log-Normal",
                                "Cauchy",  "half-Cauchy",
                                "Student t", "half-Student t", "log-Student t",
                                "Exponential", "flat (implicit)")),
    !!! prior_standat(prior_class_mean, "prior_class_mean",
                      valid = c("Normal", "Cauchy", "Student t", "flat (implicit)")),
    !!! prior_standat(prior_class_sd, "prior_class_sd",
                      valid = c("Normal", "half-Normal", "log-Normal",
                                "Cauchy",  "half-Cauchy",
                                "Student t", "half-Student t", "log-Student t",
                                "Exponential", "flat (implicit)")),
    prior_het_type = switch(prior_het_type,
                            sd = 1, var = 2, prec = 3)
    )

  # Standard pars to monitor
  pars <- c("mu", "beta", "d",
            "log_lik", "resdev",
            "lp__")

  if (has_ipd && !is_survival) pars <- c(pars, "fitted_ipd")
  if (has_agd_arm && !is_survival) pars <- c(pars, "fitted_agd_arm")
  if (has_agd_contrast) pars <- c(pars, "fitted_agd_contrast")

  # Monitor heterogeneity SD and study deltas if RE model
  if (trt_effects == "random") {
    pars <- c(pars, "tau", "delta")
  }

  # Monitor omega for node-splitting model
  if (consistency == "nodesplit") {
    pars <- c(pars, "omega")
  }

  # Monitor cumulative integration error if using numerical integration
  if (n_int > 1 && !is_survival && int_thin > 0) {
    if (has_agd_arm) pars <- c(pars, "theta_bar_cum_agd_arm")
    if (has_agd_contrast) pars <- c(pars, "theta_bar_cum_agd_contrast")
  }

  # Monitor class effects if class effects in use
  if (class_effects == "exchangeable") {
    pars <- c(pars, "class_mean", "class_sd")
  }

  # Set adapt_delta, but respect other control arguments if passed in ...
  stanargs <- list(...)
  if ("control" %in% names(stanargs))
    stanargs$control <- purrr::list_modify(stanargs$control, adapt_delta = adapt_delta)
  else
    stanargs$control <- list(adapt_delta = adapt_delta)

  # Global option rstan_refresh for refresh
  if (!"refresh" %in% names(stanargs) && !is.null(getOption("rstan_refresh"))) {
    refresh <- getOption("rstan_refresh")
    if (!rlang::is_integerish(refresh, n = 1, finite = TRUE)) abort("Global option `rstan_refresh` must be an integer.")
    stanargs$refresh <- refresh
  }

  # Set chain_id to make CHAIN_ID available in data block
  stanargs$chain_id <- 1L

  # Call Stan model for given likelihood

  # -- Normal likelihood
  if (likelihood == "normal") {

    # Dummy prior for IPD variance when no IPD - not used, but requested by Stan data
    if (!has_ipd) prior_aux <- half_normal(1)

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_y = if (has_ipd) ipd_y$.y else numeric(),
      agd_arm_y = if (has_agd_arm) agd_arm_y$.y else numeric(),
      agd_arm_se = if (has_agd_arm) agd_arm_y$.se else numeric(),

      # Add prior for auxiliary parameter - individual-level variance
      !!! prior_standat(prior_aux, "prior_aux",
                        valid = c("Normal", "half-Normal", "log-Normal",
                                  "Cauchy",  "half-Cauchy",
                                  "Student t", "half-Student t", "log-Student t",
                                  "Exponential", "flat (implicit)")),

      # Specify link
      link = switch(link, identity = 1, log = 2)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$normal,
                                   data = standat,
                                   pars = c(pars, "sigma"))

  # -- Bernoulli/binomial likelihood (one parameter)
  } else if (likelihood %in% c("bernoulli", "binomial")) {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$binomial_1par,
                                   data = standat,
                                   pars = pars)

  # -- Bernoulli/binomial likelihood (two parameter)
  } else if (likelihood %in% c("bernoulli2", "binomial2")) {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_n = if (has_agd_arm) agd_arm_y$.n else integer(),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$binomial_2par,
                                   data = standat,
                                   pars = if (n_int > 1 && int_thin > 0) c(pars, "theta2_bar_cum") else pars)

  # -- Poisson likelihood
  } else if (likelihood == "poisson") {

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ipd_r = if (has_ipd) ipd_y$.r else integer(),
      ipd_E = if (has_ipd) ipd_y$.E else numeric(),
      agd_arm_r = if (has_agd_arm) agd_arm_y$.r else integer(),
      agd_arm_E = if (has_agd_arm) agd_arm_y$.E else numeric(),

      # Specify link
      link = switch(link, log = 1)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$poisson,
                                   data = standat,
                                   pars = pars)

  # -- Ordered multinomial likelihood
  } else if (likelihood == "ordered") {

    if (has_ipd) {
      # Determine number of categories
      ncat <- ncol(ipd_y$.r)
      # Stan model takes IPD as an integer vector of category numbers
      ipd_r_int <- apply(ipd_y$.r == 1, 1, which)
      # Determine which categories are present
      ipd_cat <- t(apply(ipd_y$.r, 1,
                         function(x) {
                           cs <- which(!is.na(x))
                           c(cs, rep(0, ncat - length(cs)))
                         }))
      ipd_ncat <- rowSums(ipd_cat > 0)
    }

    if (has_agd_arm) {
      # Determine number of categories
      if (!has_ipd) ncat <- ncol(agd_arm_y$.r)
      # Determine which categories are present
      agd_arm_cat <- t(apply(agd_arm_y$.r, 1,
                         function(x) {
                           cs <- which(!is.na(x))
                           c(cs, rep(0, ncat - length(cs)))
                         }))
      agd_arm_ncat <- rowSums(agd_arm_cat > 0)
      # Replace missing category counts with 0 (these will drop out of the likelihood)
      agd_arm_r <- agd_arm_y$.r
      agd_arm_r[is.na(agd_arm_r)] <- 0
      agd_arm_n <- rowSums(agd_arm_y$.r, na.rm = TRUE)
    }

    if (!has_ipd && !has_agd_arm) {
      abort("No IPD or AgD (arm-based) in the network. Cannot fit ordered model to contrast data only.")
    }

    standat <- purrr::list_modify(standat,
      # Add outcomes
      ncat = ncat,

      ipd_r = if (has_ipd) ipd_r_int else integer(),
      ipd_cat = if (has_ipd) ipd_cat else matrix(0, 0, ncat),
      ipd_ncat = if (has_ipd) ipd_ncat else integer(),

      agd_arm_r = if (has_agd_arm) agd_arm_r else matrix(0, 0, ncat),
      agd_arm_n = if (has_agd_arm) agd_arm_n else integer(),
      agd_arm_cat = if (has_agd_arm) agd_arm_cat else matrix(0, 0, ncat),
      agd_arm_ncat = if (has_agd_arm) agd_arm_ncat else integer(),

      # Add prior for auxiliary parameters - latent cutoffs
      !!! prior_standat(prior_aux, "prior_aux",
                        valid = c("Normal", "half-Normal", "log-Normal",
                                  "Cauchy",  "half-Cauchy",
                                  "Student t", "half-Student t", "log-Student t",
                                  "Exponential", "flat (implicit)")),

      # Specify link
      link = switch(link, logit = 1, probit = 2, cloglog = 3)
    )

    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$ordered_multinomial,
                                   data = standat,
                                   pars = c(pars, "cc"))

  # -- Parametric survival likelihoods
  } else if (likelihood %in% setdiff(valid_lhood$survival, c("mspline", "pexp"))) {

    # Pull out Surv data
    ipd_surv <- get_Surv_data(ipd_y$.Surv)
    agd_arm_surv <- get_Surv_data(agd_arm_y$.Surv)

    # Add in dummy prior_aux for exponential model - not used, but requested by Stan data
    if (likelihood %in% c("exponential", "exponential-aft")) prior_aux <- flat()

    # Add in dummy prior_aux2 if not gengamma - not used, but requested by Stan data
    if (likelihood != "gengamma") prior_aux2 <- flat()

    # Add in dummy prior_aux_reg if no aux_regression
    if (is.null(X_aux)) prior_aux_reg <- flat()

    # Check aux IDs
    if (!(has_ipd || has_agd_arm)) {
      aux_id <- integer()
    } else {
      if (likelihood %in% c("exponential", "exponential-aft")) {
        aux_id <- aux_group <- rep_len(1, ni_ipd + ni_agd_arm)  # Not used (no aux pars)
      } else if (!rlang::is_integerish(aux_id, finite = TRUE)) {
        abort("`aux_id` must be an integer vector identifying the auxiliary parameter for each observation in IPD and AgD (arm-based) data.")
      }
    }

    # Set aux_int
    aux_int <- !is.null(X_aux) && max(aux_group) == nrow(X_aux)

    # Set flag for aux regression including main treatment effects
    aux_reg_trt <- any(grepl("^\\.trt[^:]", colnames(X_aux)))

    standat <- purrr::list_modify(standat,
                                  # AgD arm IDs
                                  agd_arm_arm = agd_arm_arm,

                                  # Auxiliary IDs for shape parameters
                                  #aux_by = length(aux_id) == ni_ipd + ni_agd_arm*n_int,
                                  aux_int = aux_int,
                                  aux_id = aux_id,
                                  aux_group = aux_group,

                                  # Aux regression
                                  nX_aux = if (!is.null(X_aux)) ncol(X_aux) else 0,
                                  X_aux = if (!is.null(X_aux))  X_aux else matrix(0, length(aux_id), 0),
                                  aux_reg_trt = aux_reg_trt,

                                  # Add outcomes
                                  ipd_time = ipd_surv$time,
                                  ipd_start_time = ipd_surv$start_time,
                                  ipd_delay_time = ipd_surv$delay_time,
                                  ipd_status = ipd_surv$status,

                                  agd_arm_time = agd_arm_surv$time,
                                  agd_arm_start_time = agd_arm_surv$start_time,
                                  agd_arm_delay_time = agd_arm_surv$delay_time,
                                  agd_arm_status = agd_arm_surv$status,

                                  # Specify survival distribution
                                  dist = switch(likelihood,
                                                exponential = 1,
                                                weibull = 2,
                                                gompertz = 3,
                                                `exponential-aft` = 4,
                                                `weibull-aft` = 5,
                                                lognormal = 6,
                                                loglogistic = 7,
                                                gamma = 8,
                                                gengamma = 9),

                                  # Specify link
                                  link = switch(link, log = 1),

                                  # Add prior for aux_regression coefs
                                  !!! prior_standat(prior_aux_reg, "prior_aux_reg",
                                                    valid = c("Normal", "half-Normal", "log-Normal",
                                                              "Cauchy",  "half-Cauchy",
                                                              "Student t", "half-Student t", "log-Student t",
                                                              "Exponential", "flat (implicit)"))
    )

    # Add in priors for auxiliary parameters
    if (likelihood != "gengamma") {
      standat <- purrr::list_modify(standat,
                                     # Specify prior on shape parameters
                                     !!! prior_standat(prior_aux, "prior_aux",
                                                       valid = c("Normal", "half-Normal", "log-Normal",
                                                                 "Cauchy",  "half-Cauchy",
                                                                 "Student t", "half-Student t", "log-Student t",
                                                                 "Exponential", "flat (implicit)")),

                                    # Dummy prior details for aux2
                                    !!! prior_standat(prior_aux2, "prior_aux2",
                                                      valid = "flat (implicit)"))
    } else {
      standat <- purrr::list_modify(standat,
                                     # Specify prior on sigma parameters
                                     !!! prior_standat(prior_aux$sigma, "prior_aux",
                                                       valid = c("Normal", "half-Normal", "log-Normal",
                                                                 "Cauchy",  "half-Cauchy",
                                                                 "Student t", "half-Student t", "log-Student t",
                                                                 "Exponential", "flat (implicit)")),

                                     # Specify prior on k parameters
                                     !!! prior_standat(prior_aux$k, "prior_aux2",
                                                       valid = c("Normal", "half-Normal", "log-Normal",
                                                                 "Cauchy",  "half-Cauchy",
                                                                 "Student t", "half-Student t", "log-Student t",
                                                                 "Exponential", "flat (implicit)")))
    }



    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$survival_param,
                                   data = standat,
                                   # Monitor auxiliary parameters
                                   pars =
                                     if (likelihood %in% c("exponential", "exponential-aft")) pars
                                     else if (likelihood == "lognormal") c(pars, "sdlog", "beta_aux", "d_aux")
                                     else if (likelihood == "gengamma") c(pars, "sigma", "k", "beta_aux", "d_aux")
                                     else c(pars, "shape", "beta_aux", "d_aux")
                                   )


  # -- Flexible parametric survival likelihoods (splines, piecewise exponential)
  } else if (likelihood %in% c("mspline", "pexp")) {

    # Pull out Surv data
    ipd_surv <- get_Surv_data(ipd_y$.Surv)
    agd_arm_surv <- get_Surv_data(agd_arm_y$.Surv)

    # Number of spline coefficients
    n_scoef <- ncol(basis[[1]])

    try_update <- function(..., study) {
      withCallingHandlers(update(...),
                          error = function(e) {
                            abort(glue::glue("Could not create spline basis for study {glue::double_quote(study)}."), parent = e)
                          },
                          warning = function(w) {
                            warn(glue::glue("In study {glue::double_quote(study)}: ", conditionMessage(w)))
                            rlang::cnd_muffle(w)
                          })
    }

    # Evaluate splines
    if (has_ipd) {
      ipd_time <- ipd_itime <- ipd_start_itime <- ipd_delay_itime <- matrix(nrow = length(ipd_surv$time), ncol = n_scoef)

      for (s in unique(ipd_s_t_all$.study)) {
        s_chr <- gsub("^\\.study", "", x_names[s])
        ipd_time[ipd_s_t_all$.study == s, ] <- try_update(basis[[s]], x = ipd_surv$time[ipd_s_t_all$.study == s], study = s_chr)
        ipd_itime[ipd_s_t_all$.study == s, ] <- try_update(basis[[s]], x = ipd_surv$time[ipd_s_t_all$.study == s], integral = TRUE, study = s_chr)
        ipd_start_itime[ipd_s_t_all$.study == s, ] <- try_update(basis[[s]], x = ipd_surv$start_time[ipd_s_t_all$.study == s], integral = TRUE, study = s_chr)
        ipd_delay_itime[ipd_s_t_all$.study == s, ] <- try_update(basis[[s]], x = ipd_surv$delay_time[ipd_s_t_all$.study == s], integral = TRUE, study = s_chr)
      }

    } else {
      ipd_time <- ipd_itime <- ipd_start_itime <- ipd_delay_itime <- matrix(nrow = 0, ncol = n_scoef)
    }

    if (has_agd_arm) {
      agd_arm_time <- agd_arm_itime <- agd_arm_start_itime <- agd_arm_delay_itime <- matrix(nrow = length(agd_arm_surv$time), ncol = n_scoef)

      for (s in unique(agd_arm_s_t_all$.study)) {
        s_chr <- gsub("^\\.study", "", x_names[s])
        agd_arm_time[agd_arm_s_t_all$.study == s, ] <- try_update(basis[[s]], x = agd_arm_surv$time[agd_arm_s_t_all$.study == s], study = s_chr)
        agd_arm_itime[agd_arm_s_t_all$.study == s, ] <- try_update(basis[[s]], x = agd_arm_surv$time[agd_arm_s_t_all$.study == s], integral = TRUE, study = s_chr)
        agd_arm_start_itime[agd_arm_s_t_all$.study == s, ] <- try_update(basis[[s]], x = agd_arm_surv$start_time[agd_arm_s_t_all$.study == s], integral = TRUE, study = s_chr)
        agd_arm_delay_itime[agd_arm_s_t_all$.study == s, ] <- try_update(basis[[s]], x = agd_arm_surv$delay_time[agd_arm_s_t_all$.study == s], integral = TRUE, study = s_chr)
      }

    } else {
      agd_arm_time <- agd_arm_itime <- agd_arm_start_itime <- agd_arm_delay_itime <- matrix(nrow = 0, ncol = n_scoef)
    }

    # Check aux IDs
    if (!(has_ipd || has_agd_arm)) {
      aux_id <- integer()
    } else {
      if (!rlang::is_integerish(aux_id, finite = TRUE))
        abort("`aux_id` must be an integer vector identifying the auxiliary parameter for each observation in IPD and AgD (arm-based) data.")
    }

    # Set aux_int
    aux_int <- !is.null(X_aux) && max(aux_group) == nrow(X_aux)

    # Get scoef prior means and weights for RW(1) prior with non-equally spaced knots
    if (!is.null(X_aux)) {
      # If aux_regression specified, single spline basis shared across network
      prior_aux_location <- list(mspline_constant_hazard(basis[[1]]))

      lscoef_weight <- list(rw1_prior_weights(basis[[1]]))

      aux_reg_trt <- any(grepl("^\\.trt[^:]", colnames(X_aux)))

    } else {
      prior_aux_location <- purrr::map(basis[unique(cbind(aux_id, study = c(ipd_s_t_all$.study, rep(agd_arm_s_t_all$.study, each = if (aux_int) n_int else 1))))[, "study"]],
                                       mspline_constant_hazard)

      lscoef_weight <- purrr::map(basis[unique(cbind(aux_id, study = c(ipd_s_t_all$.study, rep(agd_arm_s_t_all$.study, each = if (aux_int) n_int else 1))))[, "study"]],
                                  rw1_prior_weights)

      # Dummy prior for aux regression smoothing sd (not used)
      prior_aux_reg <- flat()

      aux_reg_trt <- FALSE
    }

    standat <- purrr::list_modify(standat,
                                  # AgD arm IDs
                                  agd_arm_arm = agd_arm_arm,

                                  # Aux IDs
                                  #aux_by = length(aux_id) == ni_ipd + ni_agd_arm*n_int,
                                  aux_int = aux_int,
                                  aux_id = aux_id,
                                  aux_group = aux_group,

                                  # Aux regression
                                  nX_aux = if (!is.null(X_aux)) ncol(X_aux) else 0,
                                  X_aux = if (!is.null(X_aux))  X_aux else matrix(0, length(aux_id), 0),
                                  aux_reg_trt = aux_reg_trt,

                                  # Number of spline coefficients
                                  n_scoef = n_scoef,

                                  # RW1 prior weights
                                  lscoef_weight = lscoef_weight,

                                  # Add outcomes
                                  ipd_time = ipd_time,
                                  ipd_itime = ipd_itime,
                                  ipd_start_itime = ipd_start_itime,
                                  ipd_delay_itime = ipd_delay_itime,
                                  ipd_delayed = ipd_surv$delay_time > 0,
                                  ipd_status = ipd_surv$status,

                                  agd_arm_time = agd_arm_time,
                                  agd_arm_itime = agd_arm_itime,
                                  agd_arm_start_itime = agd_arm_start_itime,
                                  agd_arm_delay_itime = agd_arm_delay_itime,
                                  agd_arm_delayed = agd_arm_surv$delay_time > 0,
                                  agd_arm_status = agd_arm_surv$status,

                                  # Specify link
                                  link = switch(link, log = 1),

                                  # Add prior for spline smoothing sd
                                  !!! prior_standat(prior_aux, "prior_hyper",
                                                    valid = c("Normal", "half-Normal", "log-Normal",
                                                              "Cauchy",  "half-Cauchy",
                                                              "Student t", "half-Student t", "log-Student t",
                                                              "Exponential", "flat (implicit)")),

                                  # Add prior for aux_regression smoothing sd
                                  !!! prior_standat(prior_aux_reg, "prior_reg_hyper",
                                                    valid = c("Normal", "half-Normal", "log-Normal",
                                                              "Cauchy",  "half-Cauchy",
                                                              "Student t", "half-Student t", "log-Student t",
                                                              "Exponential", "flat (implicit)"))
    )

    # Add prior mean for spline coefficients
    standat$prior_aux_location <- prior_aux_location

    # Monitor spline coefficients
    stanargs <- purrr::list_modify(stanargs,
                                   object = stanmodels$survival_mspline,
                                   data = standat,
                                   pars = c(pars,
                                            # "lscoef", "u_aux"
                                            "scoef", "beta_aux", "d_aux", "sigma", "sigma_beta"))

  } else {
    abort(glue::glue('"{likelihood}" likelihood not supported.'))
  }

  # Call sampling, managing warnings for integration checks if required
  if (n_int > 1 && int_check) {
    rhat_warn_a <- FALSE
    bulk_ess_warn_a <- FALSE
    tail_ess_warn_a <- FALSE
    stanfit <- withCallingHandlers(do.call(rstan::sampling, stanargs),
                                   warning = function(w) {
                                     m <- conditionMessage(w)
                                     if (grepl("The largest R-hat is", w, fixed = TRUE)) {
                                       rhat_warn_a <<- TRUE
                                       rlang::cnd_muffle(w)
                                     }
                                     if (grepl("Bulk Effective Sample(s?) Size \\(ESS\\) is too low", w)) {
                                       bulk_ess_warn_a <<- TRUE
                                       rlang::cnd_muffle(w)
                                     }
                                     if (grepl("Tail Effective Sample(s?) Size \\(ESS\\) is too low", w)) {
                                       tail_ess_warn_a <<- TRUE
                                       rlang::cnd_muffle(w)
                                     }
                                   })

    # Check Rhat, neff within chains with same n_int
    if (stanfit@mode == 0) {
      if (stanfit@sim$chains < nchains) {
        warn("Cannot check integration error, some chains failed to return samples.")
      } else if (any(rhat_warn_a, bulk_ess_warn_a, tail_ess_warn_a)) {
        # Only do these checks if Rhat/ESS warnings are thrown for all chains
        # combined, since these are only needed to disambiguate between low
        # iters and low n_int, and computation takes time

        sims <- as.array(stanfit)
        sims_nint1 <- sims[ , nint_vec == unique(n_int)[1], , drop = FALSE]
        sims_nint2 <- sims[ , nint_vec == unique(n_int)[2], , drop = FALSE]

        if (rhat_warn_a) {
          rhat_1 <- apply(sims_nint1, 3, rstan::Rhat)
          rhat_2 <- apply(sims_nint2, 3, rstan::Rhat)
          rhat_w <- pmax(rhat_1, rhat_2, na.rm = TRUE)
          rhat_warn_w <- any(rhat_w > 1.05, na.rm = TRUE)

          # Warnings within chains with same n_int are due to low iter
          if (rhat_warn_w) {
            warning("The largest R-hat is ", round(max(rhat_w), digits = 2),
                    ", indicating chains have not mixed.\n", "Running the chains for more iterations may help. See\n",
                    "https://mc-stan.org/misc/warnings.html#r-hat",
                    call. = FALSE)
          } else {
            # Otherwise warnings across all chains are low n_int
            rhat_a <- apply(sims, 3, rstan::Rhat)
            warn(paste0("The largest R-hat is ", round(max(rhat_a), digits = 2),
                        ", indicating that integration has not converged.\n",
                        "Increase the number of integration points `n_int` in add_integration()."),
                 class = "int_check_rhat")
          }
        }

        if (bulk_ess_warn_a) {
          bulk_ess_1 <- apply(sims_nint1, 3, rstan::ess_bulk) / ncol(sims_nint1)
          bulk_ess_2 <- apply(sims_nint2, 3, rstan::ess_bulk) / ncol(sims_nint2)
          bulk_ess_w <- pmin(bulk_ess_1, bulk_ess_2, na.rm = TRUE)
          bulk_ess_warn_w <- any(bulk_ess_w < 100, na.rm = TRUE)

          if (bulk_ess_warn_w) {
            warning("Bulk Effective Sample Size (ESS) is too low, ",
                    "indicating posterior means and medians may be unreliable.\n",
                    "Running the chains for more iterations may help. See\n",
                    "https://mc-stan.org/misc/warnings.html#bulk-ess",
                    call. = FALSE)
          } else {
            warn(paste0("Bulk Effective Sample Size (ESS) is too low, indicating that integration may not have converged.\n",
                        "Increase the number of integration points `n_int` in add_integration()."),
                 class = "int_check_essb")
          }
        }

        if (tail_ess_warn_a) {
          tail_ess_1 <- apply(sims_nint1, 3, rstan::ess_tail) / ncol(sims_nint1)
          tail_ess_2 <- apply(sims_nint2, 3, rstan::ess_tail) / ncol(sims_nint2)
          tail_ess_w <- pmin(tail_ess_1, tail_ess_2, na.rm = TRUE)
          tail_ess_warn_w <- any(tail_ess_w < 100, na.rm = TRUE)

          if (tail_ess_warn_w) {
            warning("Tail Effective Sample Size (ESS) is too low, indicating ",
                    "posterior variances and tail quantiles may be unreliable.\n",
                    "Running the chains for more iterations may help. See\n",
                    "https://mc-stan.org/misc/warnings.html#tail-ess",
                    call. = FALSE)
          } else {
            warn(paste0("Tail Effective Sample Size (ESS) is too low, indicating that integration may not have converged.\n",
                        "Increase the number of integration points `n_int` in add_integration()."),
                 class = "int_check_esst")
          }
        }
      }
    }

  } else {
    stanfit <- do.call(rstan::sampling, stanargs)
  }

  # Set readable parameter names in the stanfit object
  fnames_oi <- stanfit@sim$fnames_oi
  x_names_sub <- gsub("^(\\.study|\\.trt|\\.trtclass|\\.contr)", "", x_names)

  fnames_oi[grepl("^mu\\[[0-9]+\\]$", fnames_oi)] <- paste0("mu[", x_names_sub[col_study], "]")
  fnames_oi[grepl("^d\\[[0-9]+\\]$", fnames_oi)] <- paste0("d[", x_names_sub[col_trt], "]")
  fnames_oi[grepl("^beta\\[[0-9]+\\]$", fnames_oi)] <- paste0("beta[", x_names[col_reg], "]")
  if (!is.null(X_aux)) {
    if (likelihood %in% c("mspline", "pexp")) {
      fnames_oi[grepl("^beta_aux\\[", fnames_oi)] <- paste0("beta_aux[", rep(colnames(X_aux), times = n_scoef-1), ", ", rep(1:(n_scoef-1), each = ncol(X_aux)), "]")
      if (aux_reg_trt) {
        fnames_oi[grepl("^sigma_beta\\[", fnames_oi)] <- c("sigma_beta[.trt]",
                                                           if (ncol(X_aux) > n_trt) paste0("sigma_beta[", colnames(X_aux[, -(1:n_trt), drop = FALSE]), "]")
                                                           else NULL)
        fnames_oi[grepl("^d_aux\\[", fnames_oi)] <- paste0("d_aux[", rep(x_names_sub[col_trt], times = n_scoef-1), ", ", rep(1:(n_scoef-1), each = n_trt-1), "]")
      } else {
        fnames_oi[grepl("^sigma_beta\\[", fnames_oi)] <- paste0("sigma_beta[", colnames(X_aux), "]")
      }
    } else if (likelihood == "gengamma") {
      fnames_oi[grepl("^beta_aux\\[", fnames_oi)] <- paste0("beta_aux[", rep(colnames(X_aux), times = 2), ", ", rep(c("sigma", "k"), each = ncol(X_aux)), "]")
      if (aux_reg_trt) fnames_oi[grepl("^d_aux\\[", fnames_oi)] <- paste0("d_aux[", rep(x_names_sub[col_trt], times = 2), ", ", rep(c("sigma", "k"), each = n_trt-1), "]")
    } else {
      fnames_oi[grepl("^beta_aux\\[", fnames_oi)] <- paste0("beta_aux[", colnames(X_aux), "]")
      if (aux_reg_trt) fnames_oi[grepl("^d_aux\\[", fnames_oi)] <- paste0("d_aux[", x_names_sub[col_trt], "]")
    }
  }
  fnames_oi <- gsub("tau[1]", "tau", fnames_oi, fixed = TRUE)
  fnames_oi <- gsub("omega[1]", "omega", fnames_oi, fixed = TRUE)


  if (likelihood == "ordered") {
    if (has_ipd) l_cat <- colnames(ipd_y$.r)[-1]
    else if (has_agd_arm) l_cat <- colnames(agd_arm_y$.r)[-1]
    fnames_oi[grepl("^cc\\[[0-9]+\\]$", fnames_oi)] <- paste0("cc[", l_cat, "]")
  }

  stanfit@sim$fnames_oi <- fnames_oi

  return(stanfit)
}

#' Random effects structure
#'
#' Use `RE_cor` to generate the random effects correlation matrix, under the
#' assumption of common heterogeneity variance (i.e. all within-study
#' correlations are 0.5). Use `which_RE` to return a vector of IDs for the RE
#' deltas (0 means no RE delta on this arm).
#'
#' @param study A vector of study IDs (integer, character, or factor)
#' @param trt A factor vector of treatment codes (or coercible as such), with
#'   first level indicating the reference treatment
#' @param contrast A logical vector, of the same length as `study` and `trt`,
#'   indicating whether the corresponding data are in contrast rather than arm
#'   format.
#' @param type Character string, whether to generate RE structure under the
#'   "reference treatment" parameterisation, or the "baseline shift"
#'   parameterisation.
#'
#' @return For `RE_cor()`, a correlation matrix of dimension equal to the number
#'   of random effects deltas (excluding those that are set equal to zero).
#'
#'   For `which_RE()`, an integer vector of IDs indexing the rows and columns of
#'   the correlation matrix returned by `RE_cor()`.
#' @export
#' @aliases RE_cor
#' @rdname random_effects
#'
#' @examples
#' RE_cor(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)))
#' RE_cor(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)), type = "blshift")
RE_cor <- function(study, trt, contrast, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (!is.logical(contrast) || is.matrix(contrast))
    abort("`contrast` must be a logical vector.")
  if (length(study) != length(trt) || length(study) != length(contrast))
    abort("`study`, `trt`, and `contrast` must be the same length.")
  type <- rlang::arg_match(type)


  reftrt <- levels(trt)[1]
  if (type == "reftrt") {
    # Treat contrast rows as non ref trt arms (since they always have REs)
    nonref <- trt != reftrt | contrast
    nRE <- sum(nonref)  # RE for each non ref trt arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[nonref]
    trt <- trt[nonref]
  } else if (type == "blshift") {
    # Consider baseline arm to be first by treatment order
    # Treat contrast rows as non baseline arms (since they always have REs)
    nonbl <- tibble::tibble(study, trt, contrast) %>%
      dplyr::group_by(.data$study) %>%
      dplyr::mutate(nonbl = .data$contrast | .data$trt != sort(.data$trt)[1] | (duplicated(.data$trt) & .data$trt != reftrt)) %>%
      dplyr::pull(.data$nonbl)
    nRE <- sum(nonbl)  # RE for each non baseline arm
    Rho <- matrix(0, nrow = nRE, ncol = nRE)
    study <- study[nonbl]
    trt <- trt[nonbl]
  }

  diag(Rho) <- 1
  for (i in 1:(nRE - 1)) {
    for (j in (i + 1):nRE) {
      if (study[i] == study[j]) {
        Rho[i, j] <- 0.5
        Rho[j, i] <- 0.5
      }
    }
  }

  return(Rho)
}

#' @rdname random_effects
#' @aliases which_RE
#' @export
#' @examples
#' which_RE(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)))
#' which_RE(smoking$studyn, smoking$trtn, contrast = rep(FALSE, nrow(smoking)), type = "blshift")
which_RE <- function(study, trt, contrast, type = c("reftrt", "blshift")) {
  if (!is.numeric(study) && !is.character(study) && !is.factor(study) || is.matrix(study)) {
    abort("`study` must be a vector, either numeric, character, or factor.")
  }
  if (!is.factor(trt)) {
    trt <- tryCatch(as.factor(trt),
                    error = function(e) {
                      abort("`trt` must be a factor (or coercible to factor).")
                    }, finally = inform("Coerced `trt` to factor."))
  }
  if (!is.logical(contrast) || is.matrix(contrast))
    abort("`contrast` must be a logical vector.")
  if (length(study) != length(trt) || length(study) != length(contrast))
    abort("`study`, `trt`, and `contrast` must be the same length.")
  type <- rlang::arg_match(type)

  n_i <- length(study)
  id <- rep(0, n_i)

  reftrt <- levels(trt)[1]
  if (type == "reftrt") {
    trt_nonref <- trt != reftrt | contrast
    id[trt_nonref] <- 1:sum(trt_nonref)
  } else if (type == "blshift") {
    # Consider baseline arm to be first by treatment order
    non_bl <- tibble::tibble(study, trt, contrast) %>%
      dplyr::group_by(.data$study) %>%
      dplyr::mutate(non_bl = .data$contrast | .data$trt != sort(.data$trt)[1] | (duplicated(.data$trt) & .data$trt != reftrt)) %>%
      dplyr::pull(.data$non_bl)
    id[non_bl] <- 1:sum(non_bl)
  }

  return(id)
}

#' List of valid likelihoods for different outcomes
#' @noRd
valid_lhood <- list(binary = c("bernoulli", "bernoulli2"),
                    count = c("binomial", "binomial2"),
                    rate = "poisson",
                    continuous = "normal",
                    ordered = "ordered",
                    survival = c("exponential", "weibull", "gompertz",
                                 "exponential-aft", "weibull-aft",
                                 "lognormal", "loglogistic", "gamma", "gengamma",
                                 "mspline", "pexp"))


#' Create exchangeable class effects design vector
#' @param classes Network classes factor vector
#' @return A list, with elements `id` giving the design vector (0 = no class effect), and `label` giving the corresponding class labels
#' @noRd

which_CE <- function(classes, class_sd)   {

  # Class vector, without network reference treatment
  x <- classes[-1]

  # Identify sole occupancy classes
  solo_classes <- setdiff(levels(x)[table(x) == 1] , unlist(class_sd))

  # Set sole occupancy classes to NA (no class effects) and drop unused levels
  x <- droplevels(x, exclude = solo_classes)

  # Create numeric ID vector (0 = no class effect)
  id <- as.numeric(x)
  id[is.na(id)] <- 0

  # Create class labels
  label <- levels(x)

  return(list(id = id, label = label))
}


#' Check likelihood function, or provide default value
#'
#' @param x likelihood type as string
#' @param outcome outcome types as named list (ipd, agd_arm, agd_contrast)
#'
#' @noRd
check_likelihood <- function(x, outcome) {
  if (missing(outcome)) valid_lhood <- unlist(valid_lhood)
  else if (!is.na(outcome$ipd)) {
    valid_lhood <- valid_lhood[[outcome$ipd]]
    otype <- outcome$ipd
  } else if (!is.na(outcome$agd_arm)) {
    valid_lhood <- valid_lhood[[outcome$agd_arm]]
    otype <- outcome$agd_arm
  } else if (!is.na(outcome$agd_contrast)) {
    valid_lhood <- valid_lhood[[outcome$agd_contrast]]
    otype <- outcome$agd_contrast
  }
  else abort("No outcomes specified in network data.")

  # Choose default likelihood if NULL
  if (is.null(x) && missing(outcome)) abort("Please specify a suitable likelihood.")
  else if (is.null(x)) {
    x <- valid_lhood[1]

  # Check valid option if given
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_lhood) {
    abort(glue::glue("`likelihood` should be a character string specifying a valid likelihood.\n",
                     "Suitable options for {otype} outcomes are currently: ",
                     glue::glue_collapse(glue::double_quote(valid_lhood), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

#' Check link function, or provide default value
#'
#' @param x link function, as string
#' @param lik likelihood, as string
#'
#' @noRd
check_link <- function(x, lik) {
  valid_link <- list(normal = c("identity", "log"),
                     bernoulli = c("logit", "probit", "cloglog"),
                     bernoulli2 = c("logit", "probit", "cloglog"),
                     binomial = c("logit", "probit", "cloglog"),
                     binomial2 = c("logit", "probit", "cloglog"),
                     poisson = "log",
                     ordered = c("logit", "probit", "cloglog"),
                     exponential = "log",
                     weibull = "log",
                     gompertz = "log",
                     `exponential-aft` = "log",
                     `weibull-aft` = "log",
                     lognormal = "log",
                     loglogistic = "log",
                     gamma = "log",
                     gengamma = "log",
                     mspline = "log",
                     pexp = "log")[[lik]]

  if (is.null(x)) {
    x <- valid_link[1]
  } else if (!is.character(x) || length(x) > 1 || !tolower(x) %in% valid_link) {
    abort(glue::glue("`link` should be a character string specifying a valid link function.\n",
                     "Suitable options for a {lik} likelihood are currently: ",
                     glue::glue_collapse(glue::double_quote(valid_link), sep = ", ", last = " or "),
                     "."))
  }
  return(tolower(x))
}

#' Inverse link functions
#'
#' @param x Linear predictor values
#' @param link Character string specifying link function
#' @param ... Other parameters passed to link function
#'
#' @noRd
inverse_link <- function(x, link = c("identity", "log", "logit", "probit", "cloglog"), ...) {
  link <- rlang::arg_match(link)

  out <-
    if (link == "identity") x
    else if (link == "log") exp(x)
    else if (link == "logit") plogis(x, ...)
    else if (link == "probit") pnorm(x, ...)
    else if (link == "cloglog") 1 - exp(-exp(x))

  return(out)
}

#' Link functions
#'
#' @param x Linear predictor values
#' @param link Character string specifying link function
#' @param ... Other parameters passed to link function
#'
#' @noRd
link_fun <- function(x, link = c("identity", "log", "logit", "probit", "cloglog"), ...) {
  link <- rlang::arg_match(link)

  out <-
    if (link == "identity") x
    else if (link == "log") log(x)
    else if (link == "logit") qlogis(x, ...)
    else if (link == "probit") qnorm(x, ...)
    else if (link == "cloglog") log(-log(1 - x))

  return(out)
}

#' Get scale of outcome / linear predictor for reporting and plotting
#'
#' @param likelihood String, giving likelihood
#' @param link String, giving link function
#' @param measure String, specifying whether relative or absolute scale required
#' @param type String, specifying whether link scale or response scale required
#'
#' @return String giving the scale name, e.g. "log Odds Ratio"
#' @noRd
get_scale_name <- function(likelihood = c("normal", "bernoulli", "bernoulli2",
                                          "binomial", "binomial2", "poisson",
                                          "ordered",
                                          "exponential", "weibull", "gompertz",
                                          "exponential-aft", "weibull-aft",
                                          "lognormal", "loglogistic", "gamma",
                                          "gengamma",
                                          "mspline", "pexp"),
                           link = c("identity", "log", "logit", "probit", "cloglog"),
                           measure = c("relative", "absolute"),
                           type = c("link", "response",
                                    "survival", "hazard", "cumhaz", "mean", "median", "quantile", "rmst")) {

  likelihood <- rlang::arg_match(likelihood)
  link <- rlang::arg_match(link)
  measure <- rlang::arg_match(measure)
  type <- rlang::arg_match(type)

  check_link(link, likelihood)

  if (likelihood == "normal") {

    if (link == "identity") {
      if (measure == "relative") {
        out <- "Relative Effect"
      } else if (measure == "absolute") {
        out <- "Absolute Effect"
      }
    } else if (link == "log") {
      if (measure == "relative") {
        if (type == "link") out <- "Relative Effect (on log scale)"
        else out <- "Relative Effect (on response scale)"
      } else if (measure == "absolute") {
        if (type == "link") out <- "Absolute Effect (on log scale)"
        else out <- "Absolute Effect (on response scale)"
      }
    }

  } else if (likelihood %in% c("bernoulli", "bernoulli2", "binomial", "binomial2", "ordered")) {

    if (link == "logit") {
      if (measure == "relative") {
        if (type == "link") out <- "log Odds Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "log Odds"
        else out <- "Probability"
      }
    } else if (link == "probit") {
      if (measure == "relative") {
        if (type == "link") out <- "Probit Difference"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "Probit Probability"
        else out <- "Probability"
      }
    } else if (link == "cloglog") {
      if (measure == "relative") {
        if (type == "link") out <- "log Hazard Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "cloglog Probability"
        else out <- "Probability"
      }
    }

  } else if (likelihood == "poisson") {

    if (link == "log") {
      if (measure == "relative") {
        if (type == "link") out <- "log Rate Ratio"
        else out <- ""
      } else if (measure == "absolute") {
        if (type == "link") out <- "log Rate"
        else out <- "Rate"
      }
    }

  } else if (likelihood %in% c("exponential", "weibull", "gompertz", "mspline", "pexp")) {

    if (measure == "relative") {
      if (type == "link") out <- "log Hazard Ratio"
      else out <- ""
    } else if (measure == "absolute") {
      out <- switch(type,
                    survival = "Survival Probability",
                    hazard = "Hazard",
                    cumhaz = "Cumulative Hazard",
                    mean = "Mean Survival Time",
                    median = "Median Survival Time",
                    quantile = "Survival Time",
                    rmst = "Restricted Mean Survival Time",
                    link = "Linear Predictor")
    }

  } else if (likelihood %in% c("exponential-aft", "weibull-aft", "lognormal", "loglogistic", "gamma", "gengamma")) {

    if (measure == "relative") {
      if (type == "link") out <- "log Survival Time Ratio"
      else out <- ""
    } else if (measure == "absolute") {
      out <- switch(type,
                    survival = "Survival Probability",
                    hazard = "Hazard",
                    cumhaz = "Cumulative Hazard",
                    mean = "Mean Survival Time",
                    median = "Median Survival Time",
                    quantile = "Survival Time",
                    rmst = "Restricted Mean Survival Time",
                    link = "Linear Predictor")
    }

  } else {
    out <- ""
  }
  return(out)
}

#' Get outcome variables from internal nma_data
#'
#' @param x A data frame
#' @param o_type Outcome type
#'
#' @return A data frame with outcome variables selected
#' @noRd
get_outcome_variables <- function(x, o_type) {
  o_vars <- list(
    binary = ".r",
    rate = c(".r", ".E"),
    count = c(".r", ".n"),
    continuous = c(".y", ".se"),
    ordered = ".r",
    survival = ".Surv"
  )[[o_type]]

  return(
    dplyr::select(x, dplyr::one_of(intersect(o_vars, colnames(x))))
  )
}

#' Construct NMA formula
#'
#' @param regression Regression formula, or NULL
#' @param consistency Consistency/inconsistency model (character string)
#' @param classes Classes present? TRUE / FALSE
#' @param class_interactions Class interaction specification (character string)
#'
#' @return A formula
#' @noRd
make_nma_formula <- function(regression,
                             consistency = c("consistency", "nodesplit", "ume"),
                             classes,
                             class_interactions = c("common", "exchangeable", "independent")
                             ) {

  if (!is.null(regression) && !rlang::is_formula(regression)) abort("`regression` is not a formula")
  consistency <- rlang::arg_match(consistency)
  if (!rlang::is_bool(classes)) abort("`classes` should be TRUE or FALSE")
  if (classes && !is.null(regression)) class_interactions <- rlang::arg_match(class_interactions)

  if (!is.null(regression)) {

    # Set up treatment classes
    if (classes) {

      if (class_interactions == "common") {
        nma_formula <- do.call("substitute",
                               list(regression,
                                    list(.trt = quote(.trtclass))))

      } else if (class_interactions == "exchangeable") {
        abort('Exchangeable treatment class interactions (class_interactions = "exchangeable") not yet supported.')
      } else {
        nma_formula <- regression
      }

      # Remove any main effect of .trtclass, e.g. if user specified var*.trt
      # with common interactions, or used the .trtclass special
      nma_formula <- update.formula(nma_formula, ~. - .trtclass)

    } else {
      nma_formula <- regression
    }

    if (consistency == "ume") {
      nma_formula <- update.formula(nma_formula, ~.study + .contr + . -1)
    } else if (consistency == "nodesplit") {
      nma_formula <- update.formula(nma_formula, ~.study + .trt + .omega + . -1)
    } else {
      nma_formula <- update.formula(nma_formula, ~.study + .trt + . -1)
    }
  } else {
    if (consistency == "ume") {
      nma_formula <- ~-1 + .study + .contr
    } else if (consistency == "nodesplit") {
      nma_formula <- ~-1 + .study + .trt + .omega
    } else {
      nma_formula <- ~-1 + .study + .trt
    }
  }
}

#' Construct NMA design matrix
#'
#' @param nma_formula NMA formula, returned by [make_nma_formula()]
#' @param ipd,agd_arm,agd_contrast Data frames
#' @param agd_contrast_bl Logical vector identifying baseline rows for contrast
#'   data
#' @param xbar Named numeric vector of centering values, or NULL
#' @param consistency Consistency/inconsistency model (character string)
#' @param nodesplit Length 2 character vector giving comparison to node-split
#' @param classes Classes present? TRUE / FALSE
#' @param class_interactions Class interaction specification (character string)
#' @param newdata Providing newdata post-fitting? TRUE / FALSE
#'
#' @return A named list of three matrices: X_ipd, X_agd_arm, X_agd_contrast; and
#'   three vectors of offsets: offset_ipd, offset_agd_arm, offset_agd_contrast.
#' @noRd
make_nma_model_matrix <- function(nma_formula,
                                  dat_ipd = tibble::tibble(),
                                  dat_agd_arm = tibble::tibble(),
                                  dat_agd_contrast = tibble::tibble(),
                                  agd_contrast_bl = logical(),
                                  xbar = NULL,
                                  consistency = c("consistency", "nodesplit", "ume"),
                                  nodesplit = NULL,
                                  classes = FALSE,
                                  newdata = FALSE) {
  # Checks
  if (!rlang::is_formula(nma_formula)) abort("`nma_formula` is not a formula")
  stopifnot(is.data.frame(dat_ipd),
            is.data.frame(dat_agd_arm),
            is.data.frame(dat_agd_contrast))
  consistency <- rlang::arg_match(consistency)
  if (!rlang::is_bool(classes)) abort("`classes` should be TRUE or FALSE")
  if (nrow(dat_agd_contrast) && !rlang::is_logical(agd_contrast_bl, n = nrow(dat_agd_contrast)))
    abort("`agd_contrast_bl` should be a logical vector of length nrow(agd_contrast)")
  if (!is.null(xbar) && (
        !(rlang::is_double(xbar) || rlang::is_integer(xbar)) || !rlang::is_named(xbar)))
    abort("`xbar` should be a named numeric vector")

  if (!consistency %in% c("consistency", "ume", "nodesplit")) {
    abort(glue::glue("Inconsistency '{consistency}' model not yet supported."))
  }

  if (consistency == "nodesplit") {
    if (!rlang::is_vector(nodesplit, 2))
      abort("`nodesplit` should be a vector of length 2.")
  }

  .has_ipd <- if (nrow(dat_ipd)) TRUE else FALSE
  .has_agd_arm <- if (nrow(dat_agd_arm)) TRUE else FALSE
  .has_agd_contrast <- if (nrow(dat_agd_contrast)) TRUE else FALSE

  # Sanitise factors
  if (.has_ipd) {
    dat_ipd <- dplyr::mutate_at(dat_ipd,
      .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
      .funs = fct_sanitise)
  }
  if (.has_agd_arm) {
    dat_agd_arm <- dplyr::mutate_at(dat_agd_arm,
                            .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
                            .funs = fct_sanitise)
  }
  if (.has_agd_contrast) {
    dat_agd_contrast <- dplyr::mutate_at(dat_agd_contrast,
                            .vars = if (classes) c(".trt", ".study", ".trtclass") else c(".trt", ".study"),
                            .funs = fct_sanitise)
  }

  # Define contrasts for UME model
  if (consistency == "ume") {
    # For IPD and AgD (arm-based), take the first-ordered arm as baseline
    # So the contrast sign will always be positive
    if (.has_ipd || .has_agd_arm) {
      contrs_arm <- dplyr::bind_rows(dat_ipd, dat_agd_arm) %>%
        dplyr::distinct(.data$.study, .data$.trt) %>%
        dplyr::arrange(.data$.study, .data$.trt) %>%
        dplyr::group_by(.data$.study) %>%
        dplyr::mutate(.trt_b = dplyr::first(.data$.trt)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.contr = dplyr::if_else(.data$.trt == .data$.trt_b,
                                              "..ref..",
                                              paste0(.data$.trt, " vs. ", .data$.trt_b)),
                      .contr_sign = 1)
    } else {
      contrs_arm <- tibble::tibble()
    }

    # For AgD (contrast-based), take the specified baseline arm (with .y = NA)
    # Need to make sure to construct contrast the correct way around (d_12 instead of d_21)
    # and then make a note to change the sign if necessary
    if (.has_agd_contrast) {
      contrs_contr <- dat_agd_contrast %>%
        dplyr::distinct(.data$.study, .data$.trt) %>%  # In case integration data passed
        dplyr::left_join(dplyr::distinct(dat_agd_contrast[agd_contrast_bl, ], .data$.study, .data$.trt) %>%
                           dplyr::transmute(.data$.study, .trt_b = .data$.trt), by = ".study") %>%
        dplyr::distinct(.data$.study, .data$.trt, .data$.trt_b) %>%
        dplyr::mutate(.contr_sign = dplyr::if_else(as.numeric(.data$.trt) < as.numeric(.data$.trt_b), -1, 1),
                      .contr = dplyr::if_else(.data$.trt == .data$.trt_b,
                                              "..ref..",
                                              dplyr::if_else(.data$.contr_sign == 1,
                                                             paste0(.data$.trt, " vs. ", .data$.trt_b),
                                                             paste0(.data$.trt_b, " vs. ", .data$.trt))))
    } else {
      contrs_contr <- tibble::tibble()
    }

    # Make contrast info
    contrs_all <- dplyr::bind_rows(contrs_arm, contrs_contr)

    # Vector of all K(K-1)/2 possible contrast levels
    nt <- nlevels(contrs_all$.trt)
    ctr <- which(lower.tri(diag(nt)), arr.ind = TRUE)
    trt_lev <- levels(contrs_all$.trt)
    c_lev <- paste(trt_lev[ctr[, "row"]], trt_lev[ctr[, "col"]], sep = " vs. ")

    contrs_all <- dplyr::transmute(contrs_all,
      .data$.study, .data$.trt,
      .contr = forcats::fct_drop(factor(.data$.contr, levels = c("..ref..", c_lev))),
      .data$.contr_sign)

    # Join contrast info on to study data
    if (.has_ipd)
      dat_ipd <- dplyr::left_join(dat_ipd, contrs_all, by = c(".study", ".trt"))
    if (.has_agd_arm) {
      dat_agd_arm <- dplyr::left_join(dat_agd_arm, contrs_all, by = c(".study", ".trt"))
    }
    if (.has_agd_contrast) {
      dat_agd_contrast <- dplyr::left_join(dat_agd_contrast, contrs_all, by = c(".study", ".trt"))
    }
  }

  # Derive .omega indicator for node-splitting model
  if (consistency == "nodesplit") {
    if (.has_ipd) {
      dat_ipd <- dplyr::group_by(dat_ipd, .data$.study) %>%
        dplyr::mutate(.omega = all(nodesplit %in% .data$.trt) & .data$.trt == nodesplit[2]) %>%
        dplyr::ungroup()
    }
    if (.has_agd_arm) {
      dat_agd_arm <- dplyr::group_by(dat_agd_arm, .data$.study) %>%
        dplyr::mutate(.omega = all(nodesplit %in% .data$.trt) & .data$.trt == nodesplit[2]) %>%
        dplyr::ungroup()
    }
    if (.has_agd_contrast) {
      dat_agd_contrast <- dplyr::group_by(dat_agd_contrast, .data$.study) %>%
        dplyr::mutate(.omega = all(nodesplit %in% .data$.trt) & .data$.trt == nodesplit[2]) %>%
        dplyr::ungroup()
    }
  }

  # Construct design matrix all together then split out, so that same dummy
  # coding is used everywhere
  dat_all <- dplyr::bind_rows(dat_ipd, dat_agd_arm, dat_agd_contrast)

  # Check that required variables are present in each data set, and non-missing
  check_regression_data(nma_formula,
                        dat_ipd = dat_ipd,
                        dat_agd_arm = dat_agd_arm,
                        dat_agd_contrast = dat_agd_contrast,
                        newdata = newdata)

  # Center
  if (!is.null(xbar)) {
    dat_all[, names(xbar)] <-
      purrr::map2(dat_all[, names(xbar), drop = FALSE], xbar, ~.x - .y)
  }

  # Explicitly set contrasts attribute for key variables
  fvars <- all.vars(nma_formula)
  if (".trt" %in% fvars) stats::contrasts(dat_all$.trt) <- "contr.treatment"
  if (".trtclass" %in% fvars) stats::contrasts(dat_all$.trtclass) <- "contr.treatment"
  if (".contr" %in% fvars) stats::contrasts(dat_all$.contr) <- "contr.treatment"
  if (".omega" %in% fvars) stats::contrasts(dat_all$.omega) <- "contr.treatment"
  # .study handled separately next (not always a factor)

  # Drop study to factor to 1L if only one study (avoid contrasts need 2 or
  # more levels error)
  if (".study" %in% fvars && dplyr::n_distinct(dat_all$.study) == 1) {

    # Save study label to restore
    single_study_label <- unique(dat_all$.study)
    dat_all$.study_temp <- dat_all$.study
    dat_all$.study <- 1L

    # Fix up model formula with an intercept (if .study is main effect, which is not usually the case for aux_regression)
    if (".study" %in% colnames(attr(terms(nma_formula), "factors"))) nma_formula <- update.formula(nma_formula, ~. + 1)
  } else {
    single_study_label <- NULL
    if (".study" %in% fvars) stats::contrasts(dat_all$.study) <- "contr.treatment"
  }

  # Apply NMA formula to get design matrix
  X_all <- model.matrix(nma_formula, data = dat_all)
  offsets <- model.offset(model.frame(nma_formula, data = dat_all))
  has_offset <- !is.null(offsets)

  disc_names <- setdiff(names(attr(X_all, "contrasts")), c(".study", ".trt", ".trtclass", ".contr", ".omega"))

  if (!is.null(single_study_label)) {
    # Restore single study label and .study column
    colnames(X_all) <- stringr::str_replace(colnames(X_all),
                                            "^\\.study$",
                                            paste0(".study", single_study_label))
    dat_all <- dat_all %>%
      dplyr::mutate(.study = .data$.study_temp) %>%
      dplyr::select(-".study_temp")

    # Drop intercept column from design matrix
    if ("(Intercept)" %in% colnames(X_all)) X_all <- X_all[, -1, drop = FALSE]
  }

  # Remove columns for reference level of .trtclass
  if (classes) {
    ref_class <- levels(dat_all$.trtclass)[1]
    col_trtclass_ref <- grepl(paste0(".trtclass\\Q", ref_class, "\\E"),
                              colnames(X_all), perl = TRUE)
    X_all <- X_all[, !col_trtclass_ref, drop = FALSE]
  }

  # Remove columns for reference levels of discrete covariates
  if (length(disc_names)) for (xvar in disc_names) {
    # Check contrast type
    ctype <- attr(dat_all[[xvar]], "contrasts")
    if (is.null(ctype)) ctype <- getOption("contrasts")[if (is.ordered(dat_all[[xvar]])) "ordered" else "unordered"]

    # Get reference level for treatment/SAS contrasts
    if (ctype == "contr.treatment") {
      x_ref <-
        if (is.factor(dat_all[[xvar]])) levels(dat_all[[xvar]])[1]
        else if (is.logical(dat_all[[xvar]])) FALSE
        else levels(as.factor(dat_all[[xvar]]))[1]
    } else if (ctype == "contr.SAS") {
      x_ref <-
        if (is.factor(dat_all[[xvar]])) rev(levels(dat_all[[xvar]]))[1]
        else if (is.logical(dat_all[[xvar]])) FALSE
        else rev(levels(as.factor(dat_all[[xvar]])))[1]
    } else {
      x_ref <- NULL
    }

    # Remove reference level columns if present
    if (!is.null(x_ref)) {
      col_x_ref <- grepl(paste0("(`?)\\Q", xvar, "\\E(`?)\\Q", x_ref, "\\E"), colnames(X_all), perl = TRUE)
      X_all <- X_all[, !col_x_ref, drop = FALSE]
    }
  }

  # Remove columns for interactions with reference level of .trt or .trtclass
  ref_trt <- levels(dat_all$.trt)[1]
  regex_int_ref <- paste0("\\:\\.trt\\Q", ref_trt, "\\E$|^\\.trt\\Q", ref_trt, "\\E\\:")
  if (classes)
    regex_int_ref <- paste0(regex_int_ref, "|",
                            "\\:\\.trtclass\\Q", ref_class, "\\E$|^\\.trtclass\\Q", ref_class, "\\E\\:")
  col_int_ref <- grepl(regex_int_ref, colnames(X_all), perl = TRUE)
  X_all <- X_all[, !col_int_ref, drop = FALSE]

  # Remove global intercept column (present for aux_regression)
  X_all <- X_all[, colnames(X_all) != "(Intercept)", drop = FALSE]

  if (consistency == "ume") {
    # Set relevant entries to +/- 1 for direction of contrast, using .contr_sign
    contr_cols <- grepl("^\\.contr", colnames(X_all))
    X_all[, contr_cols] <- sweep(X_all[, contr_cols, drop = FALSE], MARGIN = 1,
                                 STATS = dat_all$.contr_sign, FUN = "*")
  }

  if (.has_ipd) {
    X_ipd <- X_all[1:nrow(dat_ipd), , drop = FALSE]
    offset_ipd <- if (has_offset) offsets[1:nrow(dat_ipd)] else NULL
  } else {
    X_ipd <- offset_ipd <- NULL
  }

  if (.has_agd_arm) {
    X_agd_arm <- X_all[nrow(dat_ipd) + 1:nrow(dat_agd_arm), , drop = FALSE]
    offset_agd_arm <- if (has_offset) offsets[nrow(dat_ipd) + 1:nrow(dat_agd_arm)] else NULL
  } else {
    X_agd_arm <- offset_agd_arm <- NULL
  }

  if (.has_agd_contrast) {
    X_agd_contrast_all <- X_all[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast), , drop = FALSE]
    offset_agd_contrast_all <-
      if (has_offset) {
        offsets[nrow(dat_ipd) + nrow(dat_agd_arm) + 1:nrow(dat_agd_contrast)]
      } else {
        NULL
      }

    # Difference out the baseline arms
    X_bl <- X_agd_contrast_all[agd_contrast_bl, , drop = FALSE]
    if (has_offset) offset_bl <- offset_agd_contrast_all[agd_contrast_bl]

    X_agd_contrast <- X_agd_contrast_all[!agd_contrast_bl, , drop = FALSE]
    offset_agd_contrast <-
      if (has_offset) {
        offset_agd_contrast_all[!agd_contrast_bl]
      } else {
        NULL
      }

    # Match non-baseline rows with baseline rows by study
    for (s in unique(dat_agd_contrast$.study)) {
      nonbl_id <- which(dat_agd_contrast$.study[!agd_contrast_bl] == s)
      bl_id <- which(dat_agd_contrast$.study[agd_contrast_bl] == s)

      bl_id <- rep_len(bl_id, length(nonbl_id))

      X_agd_contrast[nonbl_id, ] <- X_agd_contrast[nonbl_id, , drop = FALSE] - X_bl[bl_id, , drop = FALSE]
      if (has_offset) offset_agd_contrast[nonbl_id] <- offset_agd_contrast[nonbl_id] - offset_bl[bl_id]
    }

    # Remove columns for study baselines corresponding to contrast-based studies - not used
    s_contr <- unique(dat_agd_contrast$.study)
    bl_s_reg <- paste0("^\\.study(\\Q", paste0(s_contr, collapse = "\\E|\\Q"), "\\E)$")
    bl_cols <- grepl(bl_s_reg, colnames(X_agd_contrast), perl = TRUE)

    X_agd_contrast <- X_agd_contrast[, !bl_cols, drop = FALSE]
    if (.has_ipd) X_ipd <- X_ipd[, !bl_cols, drop = FALSE]
    if (.has_agd_arm) X_agd_arm <- X_agd_arm[, !bl_cols, drop = FALSE]
  } else {
    X_agd_contrast <- offset_agd_contrast <- NULL
  }

  return(list(X_ipd = X_ipd,
              X_agd_arm = X_agd_arm,
              X_agd_contrast = X_agd_contrast,
              offset_ipd = offset_ipd,
              offset_agd_arm = offset_agd_arm,
              offset_agd_contrast = offset_agd_contrast))
}

#' Extract columns used in model from data frame
#'
#' @param data Data frame
#' @param regression Regression formula or NULL
#' @param aux_regression Regression formula or NULL
#' @param label Label for data source or NULL, used for informative errors
#' @param keep Additional variables to keep in data
#'
#' @return Data frame with required columns
#' @noRd
get_model_data_columns <- function(data, regression = NULL, aux_regression = NULL, label = NULL, keep = NULL) {
  if (!is.null(label)) label <- paste(" in", label)
  regvars <- NULL
  auxregvars <- NULL

  if (!is.null(regression)) {
    regvars <- setdiff(all.vars(regression), c(".trt", ".trtclass", ".study", ".contr", ".omega"))
    badvars <- setdiff(regvars, colnames(data))
    if (length(badvars)) {
      abort(
        glue::glue("Regression variable{if (length(badvars) > 1) 's' else ''} ",
                   glue::glue_collapse(glue::double_quote(badvars), sep = ", ", last = " and "),
                   " not found", label, ".")
      )
    }
  }

  if (!is.null(aux_regression)) {
    auxregvars <- setdiff(all.vars(aux_regression), c(".trt", ".trtclass", ".study", ".contr", ".omega"))
    badvars <- setdiff(auxregvars, colnames(data))
    if (length(badvars)) {
      abort(
        glue::glue("Auxiliary regression variable{if (length(badvars) > 1) 's' else ''} ",
                   glue::glue_collapse(glue::double_quote(badvars), sep = ", ", last = " and "),
                   " not found", label, ".")
      )
    }
  }

  out <- dplyr::select(data, dplyr::starts_with("."), !! regvars, !! auxregvars, !! keep)

  # Work around dplyr::bind_rows() bug - convert .Surv column to bare matrix if present
  if (rlang::has_name(out, ".Surv")) out$.Surv <- as.matrix(out$.Surv)

  return(out)
}

#' Check data for regression
#'
#' Validates data for regression model. Checks that model matrix can be
#' constructed, and that there are no missing values.
#'
#' @param formula Model formula
#' @param dat_ipd,dat_agd_arm,dat_agd_contrast Data frames
#' @param newdata Providing newdata post-fitting? TRUE / FALSE
#'
#' @noRd
check_regression_data <- function(formula,
                                  dat_ipd = tibble::tibble(),
                                  dat_agd_arm = tibble::tibble(),
                                  dat_agd_contrast = tibble::tibble(),
                                  newdata = FALSE) {

  .has_ipd <- if (nrow(dat_ipd)) TRUE else FALSE
  .has_agd_arm <- if (nrow(dat_agd_arm)) TRUE else FALSE
  .has_agd_contrast <- if (nrow(dat_agd_contrast)) TRUE else FALSE

  # Check that required variables are present in each data set, and non-missing
  if (.has_ipd) {
    withCallingHandlers(
      X_ipd_frame <- model.frame(formula, dat_ipd, na.action = NULL),
      error = ~abort(paste0(if (newdata) "Failed to construct design matrix for `newdata`.\n" else "Failed to construct design matrix for IPD.\n", .)))

    X_ipd_has_na <- names(which(purrr::map_lgl(X_ipd_frame, ~any(is.na(.) | is.infinite(.)))))
  } else {
    X_ipd_has_na <- character(0)
  }

  if (.has_agd_arm) {
    withCallingHandlers(
      X_agd_arm_frame <- model.frame(formula, dat_agd_arm, na.action = NULL),
      error = ~abort(paste0(if (newdata) "Failed to construct design matrix for `newdata`.\n" else "Failed to construct design matrix for AgD (arm-based).\n", .)))

    X_agd_arm_has_na <- names(which(purrr::map_lgl(X_agd_arm_frame, ~any(is.na(.) | is.infinite(.)))))
  } else {
    X_agd_arm_has_na <- character(0)
  }

  if (.has_agd_contrast) {
    withCallingHandlers(
      X_agd_contrast_frame <- model.frame(formula, dat_agd_contrast, na.action = NULL),
      error = ~abort(paste0(if (newdata) "Failed to construct design matrix for `newdata`.\n" else "Failed to construct design matrix for AgD (contrast-based).\n", .)))

    X_agd_contrast_has_na <- names(which(purrr::map_lgl(X_agd_contrast_frame, ~any(is.na(.) | is.infinite(.)))))
  } else {
    X_agd_contrast_has_na <- character(0)
  }

  dat_has_na <- c(length(X_ipd_has_na) > 0,
                  length(X_agd_arm_has_na) > 0,
                  length(X_agd_contrast_has_na) > 0)
  if (any(dat_has_na)) {
    if (newdata) {
      abort(glue::glue("Variables with missing or infinite values in `newdata`: {paste(c(X_ipd_has_na, X_agd_arm_has_na, X_agd_contrast_has_na), collapse = ', ')}."))
    } else {
      abort(glue::glue(glue::glue_collapse(
        c("Variables with missing or infinite values in IPD: {paste(X_ipd_has_na, collapse = ', ')}.",
          "Variables with missing or infinite values in AgD (arm-based): {paste(X_agd_arm_has_na, collapse = ', ')}.",
          "Variables with missing or infinite values in AgD (contrast-based): {paste(X_agd_contrast_has_na, collapse = ', ')}."
        )[dat_has_na], sep = "\n")))
    }
  }

  invisible()
}

#' Check provided prior distributions
#'
#' @param x Input to check. Usually a `nma_prior` object.
#' @param list_names If `x` can be a named list of priors, the names we expect.
#'
#' @noRd
check_prior <- function(x, list_names) {
  arg <- rlang::caller_arg(x)
  if (missing(list_names)) {
    if (!inherits(x, "nma_prior"))
      abort(glue::glue("`{arg}` must be a prior distribution, see ?priors."),
            call = rlang::caller_env())
  } else {
    if (any(purrr::map_lgl(x, ~!inherits(., "nma_prior"))))
      abort(glue::glue("`{arg}` must be a named list of prior distributions, see ?priors.\n",
                       "Expecting named elements with priors for ",
                       glue::glue_collapse(list_names, sep = ", ", last = " and ", width = 30), "."),
            call = rlang::caller_env())

    nm <- setdiff(list_names, names(x))
    if (length(nm) > 0)
      abort(glue::glue("`{arg}` must be a named list of prior distributions.\n",
                       "Missing named elements with priors for ",
                       glue::glue_collapse(nm, sep = ", ", last = " and ", width = 30), "."),
            call = rlang::caller_env())
  }
}

#' Set prior details for Stan models
#'
#' @param x a `nma_prior` object
#' @param par character string, giving the Stan root parameter name (e.g.
#'   "prior_trt")
#' @param valid character vector, giving valid distributions
#'
#' @noRd
prior_standat <- function(x, par, valid){
  if (!inherits(x, "nma_prior")) abort("Not a `nma_prior` object.")

  dist <- x$dist

  if (!dist %in% valid)
    abort(glue::glue("Invalid `{par}`. Suitable distributions are: ",
                glue::glue_collapse(valid, sep = ", ", last = ", or ")))

  distn <- switch(dist,
                  `flat (implicit)` = 0,
                  Normal = , `half-Normal` = 1,
                  Cauchy = , `half-Cauchy` = 2,
                  `Student t` = , `half-Student t` = 3,
                  Exponential = 4,
                  `log-Normal` = 5,
                  `log-Student t` = 6,
                  Gamma = 7)

  out <- purrr::list_modify(c(x), dist = distn, fun = purrr::zap())
  # Set unnecessary (NA) parameters to zero. These will be ignored by Stan, but
  # need to pass rstan checks
  out[is.na(out)] <- 0
  names(out) <- paste0(par, "_", names(out))
  return(out)
}

#' Get covariance structure contrast-based data, using se on baseline arm
#'
#' @param x A data frame of agd_contrast data
#'
#' @return A list of covariance matrices, of length equal to the number of studies
#' @noRd
make_Sigma <- function(x) {
  return(unclass(
           by(x,
              forcats::fct_inorder(forcats::fct_drop(x$.study)),
              FUN = make_Sigma_block,
              simplify = FALSE)))
}

make_Sigma_block <- function(x) {
  s_ij <- x[is.na(x$.y), ".se", drop = TRUE]^2  # Covariances
  s_ii <- x[!is.na(x$.y), ".se", drop = TRUE]^2  # Variances
  narm <- length(s_ii)
  S <- matrix(s_ij, nrow = narm, ncol = narm)
  diag(S) <- s_ii
  return(S)
}

#' Get survival Surv data in consistent format
#'
#' @param x a Surv() object
#'
#' @return a list with elements time, start_time (for interval censored),
#'   delay_time (for delayed entry), and status
#' @noRd
get_Surv_data <- function(x) {
  if (!missing(x) && !is.null(x)) {

    surv_type <- attr(x, "type")

    status <- x[, "status"]

    if (surv_type == "right") {
      # Right censored
      time <- x[, "time"]
      start_time <- delay_time <- rep(0, length(time))

    } else if (surv_type == "left") {
      # Left censored
      time <- x[, "time"]
      start_time <- delay_time <- rep(0, length(time))

      # Make status indicator consistent with other types (i.e. left censored = 2)
      status[status == 0] <- 2

    } else if (surv_type %in% c("interval", "interval2")) {
      # Interval censored
      time <- start_time <- rep(0, length(x))

      # time1 is used unless status = 3 (interval censored)
      time[status < 3] <- x[status < 3, "time1"]
      time[status == 3] <- x[status == 3, "time2"]
      start_time[status == 3] <- x[status == 3, "time1"]

      delay_time <- rep(0, length(time))

    } else if (surv_type == "counting") {
      # Delayed entry
      time <- x[, "stop"]
      delay_time <- x[, "start"]
      start_time <- rep(0, length(time))
    }

  } else {
    time <- start_time <- delay_time <- numeric()
    status <- integer()
  }

  return(list(time = time,
              start_time = start_time,
              delay_time = delay_time,
              status = status))
}

#' Sanitise factor labels
#'
#' Remove forbidden chars from factor labels, replacing with "_"
#'
#' @param f A factor
#'
#' @noRd
fct_sanitise <- function(f) {
  forcats::fct_relabel(f, ~stringr::str_replace_all(., "[:\\[\\]]", "_"))
}

#' Check if formula only contains an offset
#'
#' @param f A one-sided formula, as used for regression argument in nma()
#' @noRd
is_only_offset <- function(f) {
  if (!rlang::is_formula(f, lhs = FALSE))
    abort("`f` should be a one-sided formula.")

  fterms <- terms(f)

  out <- !is.null(attr(fterms, "offset")) && length(attr(fterms, "factors")) == 0
  return(out)
}

#' Get EM variables from a model formula
#'
#' @param f A one-sided formula, as used for regression argument / output from make_nma_formula
#' @noRd
#' @return A character vector or character(0)
get_EM_vars <- function(f) {
  if (!rlang::is_formula(f, lhs = FALSE))
    abort("`f` should be a one-sided formula.")

  fnames <- colnames(attr(terms(f), "factor"))

  EM_regex <- "(^\\.trt(class)?\\:)|(\\:\\.trt(class)?$)"

  out <- stringr::str_subset(fnames, EM_regex)
  out <- unique(stringr::str_remove(out, EM_regex))
  return(out)
}

#' Make labels for data points
#'
#' @param study Study vector, coerced to character
#' @param trt Treatment vector, coerced to character
#' @param trt_b Baseline treatment vector if contrast-based, coerced to character
#'
#' @return Character vector of labels
#'
#' @noRd
make_data_labels <- function(study, trt, trt_b = NA) {
  tibble::tibble(study = study, trt = trt, trt_b = trt_b) %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(grp_n = dplyr::n(),
                  grp_id = 1:dplyr::n(),
                  vs_trt_b = dplyr::if_else(is.na(trt_b),
                                            NA_character_,
                                            paste0(" vs. ", trt_b)),
                  label = as.character(
                    dplyr::if_else(.data$grp_n > 1,
                                   glue::glue("{study}: {trt}{vs_trt_b}, {grp_id}", .na = ""),
                                   glue::glue("{study}: {trt}{vs_trt_b}", .na = "")))
    ) %>%
    dplyr::pull(.data$label)
}

#' Get design vector for auxiliary parameters
#' Currently only used for survival models, to allow the user to specify shapes
#' that vary by group (e.g. by arm) as well as study
#'
#' @param data Data frame
#' @param by Vector of variables to stratify by (character or bare column names)
#' @param add_study Add in .study factor, even if not specified?
#'
#' @noRd
get_aux_id <- function(data, by, add_study = TRUE) {
  grouped <- get_aux_by_data(data = data, by = by, add_study = add_study)
  dplyr::group_indices(grouped)
}

#' Get auxiliary parameter labels based on aux_by
#' @noRd
get_aux_labels <- function(data, by, add_study = TRUE) {
  groupdat <- get_aux_by_data(data = data, by = by, add_study = add_study) %>%
    dplyr::group_data() %>%
    dplyr::select(-".rows")

  if (ncol(groupdat) == 1) {
    paste(groupdat[[1]])
  } else if (rlang::has_name(groupdat, ".study")) {
    paste0(groupdat$.study, ": ", do.call(paste, c(dplyr::select(groupdat, -".study"), sep = ", ")))
  } else {
    do.call(paste, c(groupdat, sep = ", "))
  }
}

#' Get auxiliary group data defined by aux_by
#' @noRd
get_aux_by_data <- function(data, by, add_study = TRUE) {
  tryCatch(
    dplyr::ungroup(data) %>%
      dplyr::select(if (add_study) ".study" else NULL, !! by) %>%
      dplyr::group_by_all(),
    error = function(x) abort("`aux_by` must be a vector of variable names to stratify auxiliary parameters by.", parent = x)
  )
}

#' Determine whether auxiliary parameters need to be integrated over based on
#' aux_regression and aux_by
#' @noRd
aux_needs_integration <- function(aux_regression, aux_by) {
  (!is.null(aux_regression) && length(setdiff(colnames(attr(terms(aux_regression), "factor")), c(".study", ".trt", ".trtclass"))) > 0) ||
    (!is.null(aux_by) && length(setdiff(aux_by, c(".study", ".trt", ".trtclass"))) > 0)
}
