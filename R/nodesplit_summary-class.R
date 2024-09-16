#' The `nodesplit_summary` class
#'
#' The `nodesplit_summary` class contains posterior summary statistics for
#' node-splitting models, as a result of calling `summary()` on a
#' `nma_nodesplit` or `nma_nodesplit_df` object.
#'
#' @rdname nodesplit_summary-class
#' @name nodesplit_summary-class
#' @aliases nodesplit_summary
#'
#' @details Objects of class `nodesplit_summary` are tibble data frames, with one row
#'   for each node-split comparison and columns:
#'   \describe{
#'   \item{`trt1`, `trt2`}{Treatments forming the comparison}
#'   \item{`summary`}{A list column containing [nma_summary] objects with the
#'   posterior summaries and draws for each of the node-splitting parameters}
#'   \item{`p_value`}{Bayesian p-value for inconsistency}
#'   \item{`dic`}{A list column containing [nma_dic] objects, giving the model
#'   fit statistics}
#'   }
#'
#'   The parameters included in `summary` are:
#'   \describe{
#'   \item{`d_net`}{Network estimate from the corresponding consistency model,
#'   if available}
#'   \item{`d_dir`}{Direct estimate from the node-splitting model}
#'   \item{`d_ind`}{Indirect estimate from the node-splitting model}
#'   \item{`omega`}{Inconsistency factor \eqn{\omega = d_\mathrm{dir} -
#'   d_\mathrm{ind}}{\omega = d_dir - d_ind}}
#'   \item{`tau`}{Heterogeneity standard deviation from the node-splitting
#'   model, if a random effects model was fitted}
#'   \item{`tau_consistency`}{Heterogeneity standard deviation from the
#'   corresponding consistency model, if available and if a random effects model
#'   was fitted}
#'   }
#'
#'
NULL

#' Methods for `nodesplit_summary` objects
#'
#' The `as.data.frame()`, `as_tibble()`, and `as.tibble()` methods return the
#' node-splitting summaries in a data frame or tibble.
#'
#' @param x A `nodesplit_summary` object
#' @param ... Additional arguments passed on to other methods
#' @param digits Integer number of digits to display
#'
#' @rdname nodesplit_summary-methods
#'
#' @return A `data.frame` for `as.data.frame()`, a `tbl_df` for `as.tibble()`
#'   and `as_tibble()`.
#'
#'   The `print()` method returns `x` invisibly.
#'
#' @seealso [plot.nodesplit_summary()]
#'
#' @export
#'
print.nodesplit_summary <- function(x, ..., digits = 2) {
  if (!rlang::is_scalar_integerish(digits)) abort("`digits` must be a single integer")

  n_ns <- nrow(x)

  if (n_ns == 1) {
    cglue("Node-splitting model fitted for 1 comparison: ",
          as.character(x$trt2[1]),
          " vs. ",
          as.character(x$trt1[1]), ".")
  } else {
    cglue("Node-splitting models fitted for {n_ns} comparisons.")
  }

  for (i in 1:nrow(x)) {
    cglue("")
    if (n_ns > 1) {
      sec_header(glue::glue("Node-split {x$trt2[i]} vs. {x$trt1[i]}"))
      cglue("")
    }

    # Omit comparison details from d[]
    xsum <- x$summary[[i]]
    xsum$summary$parameter <- stringr::str_remove(xsum$summary$parameter, "\\[.*\\]")

    print(xsum, ...)
    cglue("")
    print(x$dic[[i]])
    cglue("")
    cglue("Bayesian p-value: {format.pval(x$p_value[[i]], digits = digits, eps = 10^-digits)}")
  }

  invisible(x)
}

#' @export
#' @rdname nodesplit_summary-methods
#' @param nest Whether to return a nested tibble, with the full [nma_summary]
#'   and [nma_dic] objects, or to unnest their summaries, default `FALSE`
as_tibble.nodesplit_summary <- function(x, ..., nest = FALSE) {
  if (!rlang::is_bool(nest)) abort("`nest` must be a single logical value TRUE/FALSE.")

  if (nest) { # Return underlying nested tibble
    NextMethod(...)
  } else { # Unnest summary results
    out <- x
    out$summary <- purrr::map(out$summary, tibble::as_tibble)
    out <- tidyr::unnest(out, cols = "summary")
    out$resdev <- purrr::map_dbl(out$dic, "resdev")
    out$pd <- purrr::map_dbl(out$dic, "pd")
    out$dic <- purrr::map_dbl(out$dic, "dic")
    class(out) <- setdiff(class(out), "nodesplit_summary")
    return(out)
  }
}

#' @export
#' @rdname nodesplit_summary-methods
as.tibble.nodesplit_summary <- function(x, ..., nest = FALSE) {
  return(tibble::as_tibble(x, ..., nest = nest))
}

#' @export
#' @rdname nodesplit_summary-methods
as.data.frame.nodesplit_summary <- function(x, ...) {
  return(as.data.frame(tibble::as_tibble(x, ..., nest = FALSE)))
}


#' Plots of node-splitting models
#'
#' Produce summary plots of node-splitting models
#'
#' @param x A `nodesplit_summary` object.
#' @param ... Additional arguments passed on to the underlying `ggdist` plot
#'   stat, see Details.
#' @param pars Character vector specifying the parameters to include in the
#'   plot, choices include `"d"` for the direct, indirect, and network estimates
#'   of relative effects, `"omega"` for the inconsistency factor, and `"tau"`
#'   for heterogeneity standard deviation in random effects models. Default is
#'   `"d"`.
#' @param stat Character string specifying the `ggdist` plot stat to use. The
#'   default `"dens_overlay"` is a special case, producing an overlaid density
#'   plot.
#' @param orientation Whether the `ggdist` geom is drawn horizontally
#'   (`"horizontal"`) or vertically (`"vertical"`), default `"horizontal"`.
#' @param ref_line Numeric vector of positions for reference lines, by default
#'   no reference lines are drawn.
#'
#' @details Plotting is handled by \link[ggplot2:ggplot2-package]{ggplot2} and
#'   the stats and geoms provided in the \link[ggdist:ggdist-package]{ggdist}
#'   package. As a result, the output is very flexible. Any plotting stats
#'   provided by `ggdist` may be used, via the argument `stat`. The default
#'   `"dens_overlay"` is a special exception, which uses
#'   \code{\link[ggplot2:geom_density]{ggplot2::geom_density()}}, to plot
#'   overlaid densities. Additional arguments in `...` are passed to the
#'   `ggdist` stat, to customise the output.
#'
#'   Alternative stats can be specified to produce different summaries. For
#'   example, specify `stat = "[half]eye"` to produce (half) eye plots, or `stat
#'   = "pointinterval"` to produce point estimates and credible intervals.
#'
#'   A full list of options and examples is found in the `ggdist` vignette
#'   `vignette("slabinterval", package = "ggdist")`.
#'
#'   A `ggplot` object is returned which can be further modified through the
#'   usual \link[ggplot2:ggplot2-package]{ggplot2} functions to add further
#'   aesthetics, geoms, themes, etc.
#'
#' @return A `ggplot` object.
#' @seealso [summary.nma_nodesplit_df()]
#'
#' @export
#'
#' @template ex_smoking_nma_re_nodesplit_example
#' @examples \donttest{
#' # Summarise the node-splitting results
#' (smk_nodesplit_summary <- summary(smk_fit_RE_nodesplit))
#'
#' # Plot the node-splitting results
#' plot(smk_nodesplit_summary)
#'
#' # Plot the inconsistency factors instead, change the plot stat to half-eye,
#' # and add a reference line at 0
#' plot(smk_nodesplit_summary, pars = "omega", stat = "halfeye", ref_line = 0)
#'
#' # Plot a comparison of the heterogeneity under the node-split models vs.
#' # the consistency model
#' plot(smk_nodesplit_summary, pars = "tau")
#' }
plot.nodesplit_summary <- function(x, ...,
                                   pars = "d",
                                   stat = "dens_overlay",
                                   orientation = c("horizontal", "vertical", "y", "x"),
                                   ref_line = NA_real_) {
  # Checks
  if (!is.character(pars))
    abort("`pars` must be a character vector of parameters to include.")

  if (!rlang::is_string(stat))
    abort("`stat` should be a character string specifying the name of a ggdist stat.")

  if (stat == "dens_overlay") {
    # Special case, for the classic overlaid density plot
    tb_geom <- getExportedValue("ggdist", "stat_slab")
  } else {
    stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

    tb_geom <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
                        error = function(err) {
                          abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                                      err, sep = "\n"))
                        })
  }


  if (!is.numeric(ref_line) || !is.null(dim(ref_line)))
    abort("`ref_line` should be a numeric vector.")

  orientation <- rlang::arg_match(orientation)
  if (orientation == "x") orientation <- "vertical"
  else if (orientation == "y") orientation <- "horizontal"

  # Is a horizontal geom specified?
  horizontal <- orientation == "horizontal"

  # Get draws
  draws <- tibble::as_tibble(x, nest = TRUE)
  draws$summary <- purrr::map(draws$summary, ~tibble::as_tibble(as.matrix(.)))

  draws$summary <- purrr::map(draws$summary,
                              ~tidyr::pivot_longer(., cols = dplyr::everything(),
                                 names_to = "parameter", values_to = "value"))

  draws <- dplyr::mutate(draws, comparison = paste0(.data$trt2, " vs. ", .data$trt1)) %>%
    dplyr::select("comparison", "summary") %>%
    tidyr::unnest(cols = "summary") %>%
    dplyr::mutate(parameter = stringr::str_remove(.data$parameter, "\\[.*\\]"))

  # Filter out selected pars
  # Allow selection with d instead of d_ind, d_dir, d_net
  allpars <- c(unique(draws$parameter), "d")
  badpars <- setdiff(pars, allpars)

  if (length(badpars))
    abort(glue::glue("No parameter{if (length(badpars) > 1) 's' else ''} ",
                     glue::glue_collapse(glue::double_quote(badpars), sep = ", ", last = " or "), "."))

  # Expand d and tau to include direct/indirect/network
  if ("d" %in% pars) pars <- c(setdiff(pars, "d"), "d_dir", "d_ind", "d_net")
  if ("tau" %in% pars) pars <- c(pars, "tau_consistency")

  draws <- dplyr::filter(draws, .data$parameter %in% pars)

  draws$parameter <- nfactor(draws$parameter)

  # Construct plot
  if (horizontal) {

    p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_vline(xintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::xlab("Value")

    if (stat == "dens_overlay") {
      p <- p + ggplot2::aes(colour = .data$parameter, fill = .data$parameter) +
        ggplot2::ylab("Density")
    } else {
      p <- p + ggplot2::aes(y = .data$parameter) + ggplot2::ylab("Parameter")
    }

  } else {

    p <- ggplot2::ggplot(draws, ggplot2::aes(y = .data$value)) +
      ggplot2::geom_hline(yintercept = ref_line, na.rm = TRUE, colour = "grey60") +
      ggplot2::ylab("Value")

    if (stat == "dens_overlay") {
      p <- p + ggplot2::aes(colour = .data$parameter, fill = .data$parameter) +
        ggplot2::xlab("Density")
    } else {
      p <- p + ggplot2::aes(x = .data$parameter) + ggplot2::xlab("Parameter")
    }

  }

  if (stat == "dens_overlay") {
    p <- p + do.call(ggplot2::geom_density,
                     args = rlang::dots_list(orientation = if (horizontal) "x" else "y",
                                             ...,
                                             trim = TRUE,
                                             alpha = 0.25,
                                             .homonyms = "first")) +
      ggplot2::facet_wrap(~comparison, scales = "free")
  } else {
    p <- p + do.call(tb_geom, args = list(orientation = orientation, ...)) +
      ggplot2::facet_wrap(~comparison)
  }

  p <- p +
    theme_multinma()

  return(p)
}
