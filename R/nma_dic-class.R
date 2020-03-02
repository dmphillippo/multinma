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
#'   stat? Default `TRUE`.
#' @param stat Character string specifying the `tidybayes` plot stat to use if
#'   `show_uncertainty = TRUE`, default `"pointinterval"`.
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

  if (has_y) { # Produce dev-dev plot

    # Check resdev[] names match
    if (isFALSE(all.equal(dimnames(x$resdev_array)[3], dimnames(y$resdev_array)[3])))
      abort("Data points in `x` and `y` do not match")

    # Get y resdev samples from resdev_array
    y_resdev_post <- as.matrix.nma_summary(y$resdev_array) %>%
      tibble::as_tibble()

    if (packageVersion("tidyr") >= "1.0.0") {
      y_resdev_post <- tidyr::pivot_longer(y_resdev_post, cols = dplyr::everything(),
                                           names_to = "parameter", values_to = "resdev")
    } else {
      y_resdev_post <- tidyr::gather(key = "parameter",
                                   value = "resdev",
                                   dplyr::everything())
    }

    y_resdev_post$.label <- forcats::fct_inorder(factor(
      stringr::str_extract(y_resdev_post$parameter, "(?<=\\[).+(?=\\]$)")))

    y_resdev_post$Type <- rep(Type, each = prod(dim(x$resdev_array)[1:2]))

    if (!show_uncertainty) {
      y_resdev_post <- dplyr::group_by(y_resdev_post, .data$parameter,
                                       .data$.label, .data$Type) %>%
        dplyr::summarise(resdev_y = mean(.data$resdev))

      resdev_post <- dplyr::rename(resdev_post, resdev_x = .data$resdev)

      xy_resdev_post <- dplyr::left_join(resdev_post, y_resdev_post, by = c("parameter", ".label"))
    } else {

      if (!rlang::is_string(stat))
        abort("`stat` should be a character string specifying the name of a tidybayes stat")

      stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")
      stat <- stringr::str_remove(stat, "h$")

      stath <- paste0(stat, "h")
      tb_geomh <- tryCatch(getExportedValue("tidybayes", paste0("geom_", stath)),
                           error = function(err) {
                             abort(paste("`stat` should be a character string specifying the name of a tidybayes geom:",
                                         err, sep = "\n"))
                           })

      tb_geom <- tryCatch(getExportedValue("tidybayes", paste0("geom_", stat)),
                          error = function(err) {
                            abort(paste("`stat` should be a character string specifying the name of a tidybayes geom:",
                                        err, sep = "\n"))
                          })

      # Parse ... for arguments to point_interval()
      dots <- list(...)
      is_int_arg <- names(dots) %in% c(".width", ".point", ".interval")
      int_dots <- dots[is_int_arg]
      geom_dots <- dots[!is_int_arg]

      # Summarise resdev_post and y_resdev_post with tidybayes::point_interval()
      resdev_post <- dplyr::group_by(resdev_post, .data$parameter, .data$.label, .data$Type)
      resdev_post <- do.call(tidybayes::point_intervalh,
                             args = rlang::dots_list(.data = resdev_post,
                                                     rlang::quo(.data$resdev),
                                                     .width = c(0.66, 0.95),
                                                     .point = mean,
                                                     !!! int_dots)) %>%
        dplyr::rename(resdev_x = .data$resdev,
                      x_lower = .data$.lower,
                      x_upper = .data$.upper)

      y_resdev_post <- dplyr::group_by(y_resdev_post, .data$parameter, .data$.label, .data$Type)
      y_resdev_post <- do.call(tidybayes::point_interval,
                               args = rlang::dots_list(.data = y_resdev_post,
                                                       rlang::quo(.data$resdev),
                                                       .width = c(0.66, 0.95),
                                                       .point = mean,
                                                       !!! int_dots)) %>%
        dplyr::rename(resdev_y = .data$resdev,
                      y_lower = .data$.lower,
                      y_upper = .data$.upper)

      xy_resdev_post <- dplyr::left_join(resdev_post, y_resdev_post,
                                         by = c("parameter", ".label", "Type",
                                                ".width", ".point", ".interval"))

    }

    ulim <- max(xy_resdev_post$x_upper, xy_resdev_post$y_upper)

    p <- ggplot2::ggplot(xy_resdev_post,
                         ggplot2::aes(y = .data$resdev_y,
                                      x = .data$resdev_x,
                                      xmin = .data$x_lower, xmax = .data$x_upper,
                                      ymin = .data$y_lower, ymax = .data$y_upper,
                                      colour = .data$Type,
                                      label = .data$.label)) +
      ggplot2::labs(x = "Residual Deviance (model 1)",
                    y = "Residual Deviance (model 2)") +
      ggplot2::coord_fixed(xlim = c(0, ulim), ylim = c(0, ulim)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60")

    if (dplyr::n_distinct(xy_resdev_post$Type) > 1) {
      p <- p + ggplot2::scale_colour_viridis_d("")
    } else {
      p <- p + ggplot2::guides(colour = "none")
    }

    if (show_uncertainty) {
      p <- p +
        do.call(tb_geom, args = geom_dots) +
        do.call(tb_geomh, args = geom_dots)
    } else {
      p <- p + ggplot2::geom_point(...)
    }

    p <- p + theme_multinma()

  } else { # Produce resdev plot

    if (show_uncertainty) {
      if (!rlang::is_string(stat))
        abort("`stat` should be a character string specifying the name of a tidybayes stat")

      stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

      tb_stat <- tryCatch(getExportedValue("tidybayes", paste0("stat_", stat)),
                          error = function(err) {
                            abort(paste("`stat` should be a character string specifying the name of a tidybayes stat:",
                                        err, sep = "\n"))
                          })

      # Is a horizontal geom specified?
      horizontal <- stringr::str_ends(stat, "h")
    } else {
      horizontal <- FALSE
    }

    if (horizontal) {
      resdev_post$.label <- forcats::fct_rev(resdev_post$.label)

      p <- ggplot2::ggplot(resdev_post,
                           ggplot2::aes(y = .data$.label,
                                        x = .data$resdev)) +
        ggplot2::geom_vline(xintercept = 1, colour = "grey60") +
        ggplot2::labs(x = "Residual Deviance", y = "Data Point")

      if (dplyr::n_distinct(resdev_post$Type) > 1)
        p <- p + ggplot2::facet_grid(Type~., space = "free")

    } else {
      p <- ggplot2::ggplot(resdev_post,
                           ggplot2::aes(x = .data$.label,
                                        y = .data$resdev)) +
        ggplot2::geom_hline(yintercept = 1, colour = "grey60") +
        ggplot2::labs(y = "Residual Deviance", x = "Data Point")

      if (dplyr::n_distinct(resdev_post$Type) > 1)
        p <- p + ggplot2::facet_grid(.~Type, space = "free")
    }

    if (show_uncertainty) {
      p <- p + do.call(tb_stat,
                       args = rlang::dots_list(point_interval = tidybayes::mean_qi,
                                               ...,
                                               .homonyms = "last"))
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
