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
#' ggplot(mapping = aes(colour = Treatment)) +
#'   geom_km(ndmm_net,
#'           curve_args = list(linewidth = 0.5),
#'           cens_args = list(size = 3, shape = 124)) +
#'   facet_wrap(vars(Study)) +
#'   labs(xlab = "Time", ylab = "Survival Probability") +
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
geom_km <- function(network, ..., curve_args = list(), cens_args = list()) {
  if (!inherits(network, "nma_data"))
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")

  if (!(identical(network$outcome$ipd, "survival") || identical(network$outcome$agd_arm, "survival")))
    abort("`network` does not contain survival outcomes.")

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
      )[c("time", "n.censor", "surv", "std.err", "upper", "lower")])) %>%
    # Add S(0) = 1
    dplyr::group_modify(~dplyr::add_row(., time = 0, n.censor = 0, surv = 1, std.err = 0, upper = 1, lower = 1, .before = 0)) %>%
    dplyr::mutate(Treatment = .data$.trt, Study = .data$.study)

  # Set geom args
  curve_args <- rlang::dots_list(ggplot2::aes(x = .data$time, y = .data$surv, colour = .data$Treatment, group = interaction(.data$Study, .data$Treatment)),
                                 data = kmdat,
                                 !!! curve_args,
                                 linewidth = 0.25,
                                 .homonyms = "first")

  cens_args <- rlang::dots_list(ggplot2::aes(x = .data$time, y = .data$surv, colour = .data$Treatment, group = interaction(.data$Study, .data$Treatment)),
                                data = dplyr::filter(kmdat, .data$n.censor >= 1),
                                !!! cens_args,
                                stroke = 0.25, shape = 3,
                                .homonyms = "first")

  # Output ggplot geoms
  list(do.call(ggplot2::geom_step, args = curve_args), do.call(ggplot2::geom_point, args = cens_args))
}
