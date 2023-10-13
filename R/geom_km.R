#' Kaplan-Meier curves of survival data
#'
#' This helper function constructs a ggplot2 geom to plot Kaplan-Meier curves
#' from a network containing survival or time-to-event outcomes. This is useful
#' for overlaying the "raw" survival data on the estimated survival functions
#' created with plotted with [plot.surv_nma_summary()], but can also be used
#' standalone to plot Kaplan-Meier curves before fitting a model.
#'
#' @param network A `nma_data` network object containing survival outcomes
#' @param ... Additional arguments passed to [survival::survfit()]
#' @param transform Character string giving the transformation to apply to the
#'   KM curves before plotting. The default is `"identity"` for no
#'   transformation; other options are `"cloglog"` for \eqn{\log(-\log(S))},
#'   `"log"` for \eqn{\log(S)}, or `"cumhaz"` for the cumulative hazard
#'   \eqn{-\log(S)}.
#' @param curve_args Optional list of arguments to customise the curves plotted
#' with [ggplot2::geom_step()]
#' @param cens_args Optional list of arguments to customise the censoring marks
#' plotted with [ggplot2::geom_point()]
#'
#' @return A ggplot2 geom list that can be added to a ggplot2 plot object
#' @export
#'
#' @template ex_ndmm_network
#' @examples
#' # Plot KM curves using ggplot2
#' library(ggplot2)
#'
#' # We need to create an empty ggplot object to add the curves to
#' ggplot() + geom_km(ndmm_net)
#'
#' # Adding plotting options, facets, axis labels, and a plot theme
#' ggplot() +
#'   geom_km(ndmm_net,
#'           curve_args = list(linewidth = 0.5),
#'           cens_args = list(size = 3, shape = 124)) +
#'   facet_wrap(vars(Study)) +
#'   labs(xlab = "Time", ylab = "Survival Probability") +
#'   theme_multinma()
#'
#' # Using the transform argument to produce log-log plots (e.g. to assess the
#' # proportional hazards assumption)
#' ggplot() +
#'   geom_km(ndmm_net, transform = "cloglog") +
#'   facet_wrap(vars(Study)) +
#'   theme_multinma()
#'
#' # Using the transform argument to produce cumulative hazard plots
#' ggplot() +
#'   geom_km(ndmm_net, transform = "cumhaz") +
#'   facet_wrap(vars(Study)) +
#'   theme_multinma()
#'
#' # This function can also be used to add KM data to plots of estimated survival
#' # curves from a fitted model, in a similar manner
#' @template ex_ndmm_example
#' @examples
#' # Plot estimated survival curves, and overlay the KM data
#' \donttest{
#' plot(predict(ndmm_fit, type = "survival")) + geom_km(ndmm_net)
#' }
geom_km <- function(network, ...,
                    transform = c("identity", "cloglog", "log", "cumhaz"),
                    curve_args = list(), cens_args = list()) {
  if (!inherits(network, "nma_data"))
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")

  if (!(identical(network$outcome$ipd, "survival") || identical(network$outcome$agd_arm, "survival")))
    abort("`network` does not contain survival outcomes.")

  transform <- rlang::arg_match(transform)

  if (!is.list(curve_args))
    abort("`curve_args` must be a list of arguments to pass to `geom_step()`.")

  if (!is.list(cens_args))
    abort("`cens_args` must be a list of arguments to pass to `geom_point()`.")

  dots <- list(...)

  # Get KM fits
  kmdat <- dplyr::bind_rows(network$ipd,
                            if (has_agd_arm(network)) tidyr::unnest(network$agd_arm, cols = ".Surv") else NULL) %>%
    dplyr::group_by(.data$.study, .data$.trt) %>%
    dplyr::group_modify(~dplyr::as_tibble(unclass(
      do.call(survival::survfit, rlang::dots_list(formula = .Surv ~ 1, !!! dots,  data = ., .homonyms = "last"))
      )[c("time", "n.censor", "surv", "std.err", "upper", "lower", "cumhaz")]))

  # Add S(0) = 1
  if (transform != "cloglog") kmdat <- dplyr::group_modify(kmdat, ~dplyr::add_row(., time = 0, n.censor = 0, surv = 1, std.err = 0, upper = 1, lower = 1, cumhaz = 0, .before = 0))

  kmdat <- dplyr::mutate(kmdat, Treatment = .data$.trt, Study = .data$.study)

  aes <- purrr::list_modify(
    ggplot2::aes(x = .data$time, y = .data$surv,
                 colour = .data$Treatment,
                 group = interaction(.data$Study, .data$Treatment)),
    !!! switch(transform,
           "cloglog" = ggplot2::aes(x = log(.data$time), y = log(-log(.data$surv))),
           "log" = ggplot2::aes(y = log(.data$surv)),
           "cumhaz" = ggplot2::aes(y = .data$cumhaz))
  )

  # Set geom args
  curve_args <- rlang::dots_list(mapping = aes,
                                 data = kmdat,
                                 !!! curve_args,
                                 linewidth = 0.25,
                                 .homonyms = "first")

  cens_args <- rlang::dots_list(mapping = aes,
                                data = dplyr::filter(kmdat, .data$n.censor >= 1),
                                !!! cens_args,
                                stroke = 0.25, shape = 3,
                                .homonyms = "first")

  # Output ggplot geoms
  list(do.call(ggplot2::geom_step, args = curve_args), do.call(ggplot2::geom_point, args = cens_args))
}
