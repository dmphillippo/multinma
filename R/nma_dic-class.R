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
                         show_uncertainty = TRUE,
                         stat = "pointinterval") {
  # Checks
  has_y <- !missing(y)

  if (has_y && !inherits(y, "nma_dic"))
    abort("Second argument `y` should be a `nma_dic` object produced by dic(), or missing.")

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

    # Check resdev[] names match

  } else { # Produce resdev plot

    # Get resdev samples from resdev_array
    resdev_post <- as.matrix.nma_summary(x$resdev_array) %>%
      tibble::as_tibble()

    if (packageVersion("tidyr") >= "1.0.0") {
      resdev_post <- tidyr::pivot_longer(resdev_post, cols = dplyr::everything(),
                                   names_to = "parameter", values_to = "resdev")
    } else {
      resdev_post <- tidyr::gather(key = "parameter",
                             value = "resdev",
                             dplyr::everything())
    }

    resdev_post$.label <- forcats::fct_inorder(factor(
      stringr::str_extract(resdev_post$parameter, "(?<=\\[).+(?=\\]$)")))

    Type <- c(rep("IPD", NROW(x$pointwise$ipd)),
              rep("AgD (arm-based)", NROW(x$pointwise$agd_arm)),
              rep("AgD (contrast-based)", NROW(x$pointwise$agd_contrast)))

    resdev_post$Type <- rep(Type, each = prod(dim(x$resdev_array)[1:2]))

    if (!show_uncertainty) {
      resdev_post <- dplyr::group_by(resdev_post, .data$parameter,
                                     .data$.label, .data$Type) %>%
        dplyr::summarise(resdev = mean(.data$resdev))
    }

    if (horizontal) {
      resdev_post$.label <- forcats::fct_rev(resdev_post$.label)

      p <- ggplot2::ggplot(resdev_post,
                           ggplot2::aes(y = .data$.label,
                                        x = .data$resdev)) +
        ggplot2::geom_vline(xintercept = 1, colour = "grey60") +
        ggplot2::facet_grid(Type~., space = "free") +
        ggplot2::labs(x = "Residual Deviance", y = "Data Point")
    } else {
      p <- ggplot2::ggplot(resdev_post,
                           ggplot2::aes(x = .data$.label,
                                        y = .data$resdev)) +
        ggplot2::geom_hline(yintercept = 1, colour = "grey60") +
        ggplot2::facet_grid(.~Type, space = "free") +
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
  }

  return(p)
}
