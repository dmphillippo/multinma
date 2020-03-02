#' The nma_dic class
#'
#' The `nma_dic` class contains details of the Deviance Information Criterion
#' (DIC).
#'
#' @rdname nma_dic-class
#' @name nma_dic-class
#' @aliases nma_dic
#'
#' @details Objects of class `nma_dic` have the following components:
#'   \describe{
#'   \item{`dic`}{The DIC value}
#'   \item{`pd`}{The effective number of parameters}
#'   \item{`resdev`}{The total residual deviance}
#'   \item{`pointwise`}{A list of data frames containing the pointwise
#'   contributions for the IPD and AgD.}
#'   \item{`resdev_array`}{A 3D MCMC array \[Iterations, Chains, Parameters\] of
#'   posterior residual deviance samples.}
#'   }
#'
NULL

#' Print DIC details
#'
#' @param x An object of class [nma_dic]
#' @param digits An integer passed to [round()]
#' @param ... Ignored
#'
#' @return
#' @export
#'
#' @examples
print.nma_dic <- function(x, digits = 1, ...) {
  if (!is.numeric(digits) ||
      length(digits) > 1 ||
      trunc(digits) != digits) abort("`digits` must be a single integer.")

  n <- sum(nrow(x$pointwise$ipd),
           nrow(x$pointwise$agd_arm),
           x$pointwise$agd_contrast$n_contrast)

  cglue("Residual deviance: {round(x$resdev, digits)}", subtle(" (on {n} data points)", sep = ""))
  cglue("               pD: {round(x$pd, digits)}")
  cglue("              DIC: {round(x$dic, digits)}")
}


#' Plots of model fit diagnostics
#'
#' The `plot()` method for [nma_dic] objected produced by [dic()] produces
#' several useful diagnostic plots for checking model fit and model comparison.
#' Further detail on these plots and their interpretation is given by
#' \insertCite{TSD2;textual}{multinma}.
#'
#' @param x A `nma_dic` object
#' @param y (Optional) A second `nma_dic` object, to produce "dev-dev" plots for
#'   model comparison.
#' @param ... Additional arguments passed on to other methods
#' @param type Type of plot to produce, either `"resdev"` (the default) for
#'   plots of residual deviance contributions, or `"leverage"` for plots of
#'   leverage (pD) against residual deviance. Only `"resdev"` is supported if
#'   `y` is provided.
#' @param show_uncertainty Logical, show uncertainty with a `tidybayes` plot
#'   stat? Default `TRUE`. Only used when `type = "resdev"`.
#' @param stat Character string specifying the `tidybayes` plot stat to use if
#'   `show_uncertainty = TRUE`, default `"pointinterval"`. Only used when `type
#'   = "resdev"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
plot.nma_dic <- function(x, y, ...,
                         type = c("resdev", "pd"),
                         show_uncertainty = TRUE,
                         stat = "pointinterval") {
  # Checks
  has_y <- !missing(y)

  if (has_y && !inherits(y, "nma_dic"))
    abort("Second argument `y` should be a `nma_dic` object produced by dic(), or missing.")

  if (!has_y)
    type <- rlang::arg_match(type)

  if (!rlang::is_bool(show_uncertainty))
    abort("`show_uncertainty` should be TRUE or FALSE.")

  if (show_uncertainty) {
    if (!rlang::is_string(stat))
      abort("`stat` should be a character string specifying the name of a tidybayes stat")

    stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

    tb_geom <- tryCatch(getExportedValue("tidybayes", paste0("stat_", stat)),
                        error = function(err) {
                          abort(paste("`stat` should be a character string specifying the name of a tidybayes stat:",
                                      err, sep = "\n"))
                        })

    # Is a horizontal geom specified?
    horizontal <- stringr::str_ends(stat, "h")
  } else {
    horizontal <- FALSE
  }

  if (has_y) { # Produce dev-dev plot

  } else { # Produce resdev or leverage plots
    if (!is.null(x$pointwise$ipd)) {
      dat_ipd <- x$pointwise$ipd %>%
        dplyr::mutate(.id = 1:dplyr::n(),
                      .label = as.character(glue::glue("{.study}: {.trt}, {.id}")),
                      Type = "IPD")
    } else {
      dat_ipd <- tibble::tibble()
    }

    if (!is.null(x$pointwise$agd_arm)) {
      dat_agd_arm <- x$pointwise$agd_arm %>%
        dplyr::mutate(.label = as.character(glue::glue("{.study}: {.trt}")),
                      Type = "AgD (arm-based)")
    } else {
      dat_agd_arm <- tibble::tibble()
    }

    if (!is.null(x$pointwise$agd_contrast)) {
      dat_agd_contrast <- x$pointwise$agd_contrast %>%
        dplyr::mutate(.label = as.character(glue::glue("{.study} ({n_contrast + 1} arms)")),
                      Type = "AgD (contrast-based)")
    } else {
      dat_agd_contrast <- tibble::tibble()
    }

    dat_all <- dplyr::bind_rows(dat_ipd, dat_agd_arm, dat_agd_contrast)
    dat_all$.label <- forcats::fct_inorder(factor(dat_all$.label))

    if (type == "resdev") {
      if (horizontal) {
        dat_all$.label <- forcats::fct_rev(dat_all$.label)

        p <- ggplot2::ggplot(dat_all,
                             ggplot2::aes(y = .data$.label,
                                          x = .data$resdev)) +
          ggplot2::geom_vline(xintercept = 1, colour = "grey60") +
          ggplot2::facet_grid(Type~., scales = "free", space = "free") +
          ggplot2::labs(x = "Residual Deviance", y = "Data Point")
      } else {
        p <- ggplot2::ggplot(dat_all,
                             ggplot2::aes(x = .data$.label,
                                          y = .data$resdev)) +
          ggplot2::geom_hline(yintercept = 1, colour = "grey60") +
          ggplot2::facet_grid(.~Type, scales = "free", space = "free") +
          ggplot2::labs(y = "Residual Deviance", x = "Data Point")
      }

      if (show_uncertainty) {
        p <- p + do.call(tb_geom, args = list(...))
      } else {
        p <- p + ggplot2::geom_point(...)
      }

      p <- p + theme_multinma()

      if (!horizontal) {
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0))
      } else {
        p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0))
      }

    } else if (type == "leverage") {

    }
  }

  return(p)
}
