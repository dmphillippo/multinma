#' The nma_dic class
#'
#' The `nma_dic` class contains details of the Deviance Information Criterion
#' (DIC), produced using the [dic()] function.
#'
#' @rdname nma_dic-class
#' @name nma_dic-class
#' @aliases nma_dic
#'
#' @details Objects of class `nma_dic` have the following components:
#'   \describe{
#'   \item{`dic`}{The DIC value}
#'   \item{`pd`, `pv`}{The effective number of parameters}
#'   \item{`resdev`}{The total residual deviance}
#'   \item{`pointwise`}{A list of data frames containing the pointwise
#'   contributions for the IPD and AgD.}
#'   \item{`resdev_array`}{A 3D MCMC array \[Iterations, Chains, Parameters\] of
#'   posterior residual deviance samples.}
#'   }
#'
#' @seealso [dic()], [print.nma_dic()], [plot.nma_dic()].
#'
NULL

#' Methods for `nma_dic` objects
#'
#' The `print()` method prints details of DIC model fit statistics, computed by
#' the  [dic()] function. The `as.data.frame()`, `as_tibble()`, and `as.tibble()` methods return the
#' pointwise contributions to the DIC and $p_D$ in a data frame or tibble. The `as.array()` and `as.matrix()`
#' methods returns a 3D MCMC array (as class [mcmc_array]) or matrix of posterior draws of the residual deviances.
#'
#'
#' @param x An object of class [nma_dic]
#' @param digits Integer number of digits to display
#' @param ... Additional arguments passed on to other methods
#'
#' @rdname nma_dic-methods
#'
#' @seealso [dic()], [plot.nma_dic()]
#'
#' @return A `data.frame` for `as.data.frame()`, a `tbl_df` for `as.tibble()`
#'   and `as_tibble()`, a `matrix` for `as.matrix()`, and an `mcmc_array` for
#'   `as.array()`.
#'
#'   The `print()` method returns `x` invisibly.
#' @export
#'
print.nma_dic <- function(x, digits = 1, ...) {
  if (!rlang::is_scalar_integerish(digits)) abort("`digits` must be a single integer.")

  penalty <- attr(x, "penalty")

  n <- sum(if (rlang::has_name(x$pointwise$ipd, "df")) x$pointwise$ipd$df else nrow(x$pointwise$ipd),
           if (rlang::has_name(x$pointwise$agd_arm, "df")) x$pointwise$agd_arm$df else nrow(x$pointwise$agd_arm),
           if (rlang::has_name(x$pointwise$agd_contrast, "df")) x$pointwise$agd_contrast$df else x$pointwise$agd_contrast$n_contrast)

  cglue("Residual deviance: {round(x$resdev, digits)}", subtle(" (on {n} data points)", sep = ""))
  if (penalty == "pD")
    cglue("               pD: {round(x$pd, digits)}")
  else if (penalty == "pV")
    cglue("               pV: {round(x$pv, digits)}")
  cglue("              DIC: {round(x$dic, digits)}")

  invisible(x)
}

#' @rdname nma_dic-methods
#' @export
as.data.frame.nma_dic <- function(x, ...) {
  return(as.data.frame(dplyr::bind_rows(purrr::map(x$pointwise, ~dplyr::select(., -"observed", -"fitted")), ...)))
}

#' @rdname nma_dic-methods
#' @method as.tibble nma_dic
#' @export
as.tibble.nma_dic <- function(x, ...) {
  return(dplyr::bind_rows(purrr::map(x$pointwise, ~dplyr::select(., -"observed", -"fitted")), ...))
}

#' @rdname nma_dic-methods
#' @export
as_tibble.nma_dic <- function(x, ...) {
  return(dplyr::bind_rows(purrr::map(x$pointwise, ~dplyr::select(., -"observed", -"fitted")), ...))
}

#' @rdname nma_dic-methods
#' @export
as.array.nma_dic <- function(x, ...) {
  out <- x$resdev_array
  class(out) <-  c("mcmc_array", "array")
  return(out)
}

#' @rdname nma_dic-methods
#' @export
as.matrix.nma_dic <- function(x, ...){
  # Follow approach in rstan:::as.matrix.stanfit
  a <- as.array(x)
  names_a <- dimnames(a)
  dim_a <- dim(a)
  dim(a) <- c(dim_a[1] * dim_a[2], dim_a[3])
  dimnames(a) <- names_a[-2]
  class(a) <- "matrix"
  return(a)
}

#' Plots of model fit diagnostics
#'
#' The `plot()` method for [nma_dic] objects produced by [dic()] produces
#' several useful diagnostic plots for checking model fit and model comparison.
#' Further detail on these plots and their interpretation is given by
#' \insertCite{TSD2;textual}{multinma}.
#'
#' @param x A `nma_dic` object
#' @param y (Optional) A second `nma_dic` object, to produce "dev-dev" plots for
#'   model comparison.
#' @param ... Additional arguments passed on to other methods
#' @param type With a single `nma_dic` object, whether to show the residual
#'   deviance contribution for each data point (`"resdev"`, the default), or a
#'   leverage plot (`"leverage"`).
#' @param show_uncertainty Logical, show uncertainty with a `ggdist` plot
#'   stat? Default `TRUE`. Ignored for `type = "leverage"`.
#' @param stat Character string specifying the `ggdist` plot stat to use if
#'   `show_uncertainty = TRUE`, default `"pointinterval"`. If `y` is provided,
#'   currently only `"pointinterval"` is supported.
#' @param orientation Whether the `ggdist` geom is drawn horizontally
#'   (`"horizontal"`) or vertically (`"vertical"`). Only used for residual
#'   deviance plots, default `"vertical"`.
#' @param dic_contours Numeric vector of DIC contribution contours to show on
#'   leverage plots, default `1:4`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @details When a single `nma_dic` object is given, the default plot (`type =
#'   "resdev"`) shows the residual deviance contribution for each data point.
#'   For a good fitting model, each data point is expected to have a residual
#'   deviance equal to its degrees of freedom (typically 1, except for multi-arm
#'   trials in contrast format or multinomial outcomes); larger values indicate
#'   data points that are fit poorly by the model. A leverage plot can be
#'   produced with `type = "leverage"`, plotting the leverages (contributions to
#'   \eqn{p_D}) against the signed square root residual deviance, both
#'   standardised by the degrees of freedom. Contours for different DIC
#'   contribution cutoffs are also shown; points which lie outside the DIC=3
#'   contour are generally identified as contributing to a model's poor fit.
#'
#'   When two `nma_dic` objects are given, a "dev-dev" plot comparing the
#'   residual deviance contributions under each model is produced. Data points
#'   with residual deviance contributions lying on the line of equality are fit
#'   equally well under either model. Data points lying below the line of
#'   equality indicate better fit under the second model (`y`); conversely, data
#'   points lying above the line of equality indicate better fit under the first
#'   model (`x`). A common use case is to compare a standard consistency model
#'   (fitted using [nma()] with `consistency = "consistency"`) with an unrelated
#'   mean effects (UME) inconsistency model (fitted using [nma()] with
#'   `consistency = "ume"`), to check for potential inconsistency.
#'
#'   See \insertCite{TSD2;textual}{multinma} for further details.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples ## Smoking cessation
#' @template ex_smoking_nma_fe_example
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Compare DIC of FE and RE models
#' (smk_dic_FE <- dic(smk_fit_FE))
#' (smk_dic_RE <- dic(smk_fit_RE))   # substantially better fit
#'
#' # Plot residual deviance contributions under RE model
#' plot(smk_dic_RE)
#'
#' # Further customisation is possible using ggplot commands
#' # For example, highlighting data points with residual deviance above a certain threshold
#' plot(smk_dic_RE) +
#'   ggplot2::aes(colour = ggplot2::after_stat(ifelse(y > 1.5, "darkorange", "black"))) +
#'   ggplot2::scale_colour_identity()
#'
#' # Or by posterior probability, for example here a central probability of 0.6
#' # corresponds to a lower tail probability of (1 - 0.6)/2 = 0.2
#' plot(smk_dic_RE, .width = c(0.6, 0.95)) +
#'   ggplot2::aes(colour = ggplot2::after_stat(ifelse(ymin > 1, "darkorange", "black"))) +
#'   ggplot2::scale_colour_identity()
#'
#' # Leverage plot for RE model
#' plot(smk_dic_RE, type = "leverage")
#'
#' # Check for inconsistency using UME model
#' }
#' @template ex_smoking_nma_re_ume_example
#' @examples \donttest{
#' # Compare DIC
#' smk_dic_RE
#' (smk_dic_RE_UME <- dic(smk_fit_RE_UME))  # no difference in fit
#'
#' # Compare residual deviance contributions with a "dev-dev" plot
#' plot(smk_dic_RE, smk_dic_RE_UME)
#'
#' # By default the dev-dev plot can be a little cluttered
#' # Hiding the credible intervals
#' plot(smk_dic_RE, smk_dic_RE_UME, show_uncertainty = FALSE)
#'
#' # Changing transparency
#' plot(smk_dic_RE, smk_dic_RE_UME, point_alpha = 0.5, interval_alpha = 0.1)
#' }
plot.nma_dic <- function(x, y, ...,
                         type = c("resdev", "leverage"),
                         show_uncertainty = TRUE,
                         stat = "pointinterval",
                         orientation = c("vertical", "horizontal", "x", "y"),
                         dic_contours = 1:4) {
  # Checks
  has_y <- !missing(y)

  if (has_y && !inherits(y, "nma_dic"))
    abort("Second argument `y` should be a `nma_dic` object produced by dic(), or missing.")

  if (!rlang::is_bool(show_uncertainty))
    abort("`show_uncertainty` should be TRUE or FALSE.")

  orientation <- rlang::arg_match(orientation)
  if (orientation == "x") orientation <- "vertical"
  else if (orientation == "y") orientation <- "horizontal"

  type <- rlang::arg_match(type)
  if (type == "leverage") {
    show_uncertainty <- FALSE
    if (!rlang::is_vector(dic_contours) || !(is.numeric(dic_contours) || all(is.na(dic_contours))))
      abort("`dic_contours` must be a numeric vector.")
    if (!is.null(x$pointwise$ipd) && survival::is.Surv(x$pointwise$ipd[["observed"]]) ||
        !is.null(x$pointwise$agd_arm) && survival::is.Surv(x$pointwise$agd_arm[["observed"]]))
      abort("Leverage plots are not currently supported for survival outcomes.")
  }

  # Get resdev samples from resdev_array
  resdev_post <- as.matrix.nma_summary(x$resdev_array) %>%
    tibble::as_tibble()

  data_labels <- stringr::str_extract(colnames(resdev_post), "(?<=\\[).+(?=\\]$)")
  colnames(resdev_post) <- data_labels

  resdev_post <- tidyr::pivot_longer(resdev_post, cols = dplyr::everything(),
                                     names_to = "parameter", values_to = "resdev")

  resdev_post$.label <- forcats::fct_inorder(factor(resdev_post$parameter))

  # Make sure rows for parameters are together - behaviour change between gather() and pivot_longer()
  resdev_post$.ord <- factor(resdev_post$parameter, levels = data_labels)
  resdev_post <- dplyr::arrange(resdev_post, .data$.ord)

  Type <- c(rep("IPD", NROW(x$pointwise$ipd)),
            rep("AgD (arm-based)", NROW(x$pointwise$agd_arm)),
            rep("AgD (contrast-based)", NROW(x$pointwise$agd_contrast)))

  resdev_post$Type <- rep(Type, each = prod(dim(x$resdev_array)[1:2]))

  # Some likelihoods have data points with >1 df, these are stored in `df` column
  has_df_ipd <- rlang::has_name(x$pointwise$ipd, "df")
  has_df_agd_arm <- rlang::has_name(x$pointwise$agd_arm, "df")
  has_df_agd_contrast <- rlang::has_name(x$pointwise$agd_contrast, "df")

  has_df <- has_df_ipd || has_df_agd_arm || has_df_agd_contrast

  df <- c(if (has_df_ipd) x$pointwise$ipd$df else rep(1, NROW(x$pointwise$ipd)),
          if (has_df_agd_arm) x$pointwise$agd_arm$df else rep(1, NROW(x$pointwise$agd_arm)),
          if (has_df_agd_contrast) x$pointwise$agd_contrast$df else rep(1, NROW(x$pointwise$agd_contrast)))

  resdev_post$df <- rep(df, each = prod(dim(x$resdev_array)[1:2]))
  resdev_post$df_label <- paste0("df = ", resdev_post$df)

  if (!show_uncertainty) {
    resdev_post <- dplyr::group_by(resdev_post, .data$parameter,
                                   .data$.label, .data$Type,
                                   .data$df, .data$df_label) %>%
      dplyr::summarise(resdev = mean(.data$resdev))
  }

  if (has_y) { # Produce dev-dev plot

    # Get y resdev samples from resdev_array
    y_resdev_post <- as.matrix.nma_summary(y$resdev_array) %>%
      tibble::as_tibble()

    y_data_labels <- stringr::str_extract(colnames(y_resdev_post), "(?<=\\[).+(?=\\]$)")
    colnames(y_resdev_post) <- y_data_labels

    # Check resdev[] names match
    if (any(data_labels != y_data_labels))
      abort("Data points in `x` and `y` do not match")

    y_resdev_post <- tidyr::pivot_longer(y_resdev_post, cols = dplyr::everything(),
                                         names_to = "parameter", values_to = "resdev")

    y_resdev_post$.label <- forcats::fct_inorder(factor(y_resdev_post$parameter))

    # Make sure rows for parameters are together - behaviour change between gather() and pivot_longer()
    y_resdev_post$.ord <- factor(y_resdev_post$parameter, levels = y_data_labels)
    y_resdev_post <- dplyr::arrange(y_resdev_post, .data$.ord)

    y_resdev_post$Type <- rep(Type, each = prod(dim(x$resdev_array)[1:2]))

    if (!show_uncertainty) {
      y_resdev_post <- dplyr::group_by(y_resdev_post, .data$parameter,
                                       .data$.label, .data$Type) %>%
        dplyr::summarise(resdev_y = mean(.data$resdev))

      resdev_post <- dplyr::rename(resdev_post, resdev_x = "resdev")

      xy_resdev_post <- dplyr::left_join(resdev_post, y_resdev_post,
                                         by = c("parameter", ".label", "Type")) %>%
        # Plot IPD points underneath for better clarity
        dplyr::arrange(dplyr::desc(.data$Type), .data$parameter, .data$.label)
    } else {

      if (!rlang::is_string(stat))
        abort("`stat` should be a character string specifying the name of a ggdist geom")

      stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

      if (stat != "pointinterval")
        warn(glue::glue('Currently only stat = "pointinterval" is supported when `y` is given.',
                        'Results may be unexpected with stat = "{stat}"!',
                        .sep = "\n"))

      tb_geom <- tryCatch(getExportedValue("ggdist", paste0("geom_", stat)),
                          error = function(err) {
                            abort(paste("`stat` should be a character string specifying the name of a ggdist geom:",
                                        err, sep = "\n"))
                          })

      # Parse ... for arguments to point_interval()
      dots <- list(...)
      is_int_arg <- names(dots) %in% c(".width", ".point", ".interval")
      int_dots <- dots[is_int_arg]
      geom_dots <- dots[!is_int_arg]

      # Summarise resdev_post and y_resdev_post with ggdist::point_interval()
      resdev_post <- dplyr::group_by(resdev_post, .data$parameter, .data$.label, .data$Type)
      resdev_post <- do.call(ggdist::point_interval,
                             args = rlang::dots_list(.data = resdev_post,
                                                     rlang::quo(.data$resdev),
                                                     .width = c(0.66, 0.95),
                                                     .point = mean,
                                                     !!! int_dots,
                                                     .homonyms = "last")) %>%
        dplyr::rename(resdev_x = "resdev",
                      x_lower = ".lower",
                      x_upper = ".upper")

      y_resdev_post <- dplyr::group_by(y_resdev_post, .data$parameter, .data$.label, .data$Type)
      y_resdev_post <- do.call(ggdist::point_interval,
                               args = rlang::dots_list(.data = y_resdev_post,
                                                       rlang::quo(.data$resdev),
                                                       .width = c(0.66, 0.95),
                                                       .point = mean,
                                                       !!! int_dots,
                                                       .homonyms = "last")) %>%
        dplyr::rename(resdev_y = "resdev",
                      y_lower = ".lower",
                      y_upper = ".upper")

      xy_resdev_post <- dplyr::left_join(resdev_post, y_resdev_post,
                                         by = c("parameter", ".label", "Type",
                                                ".width", ".point", ".interval"))

    }

    if (show_uncertainty) {
      ulim <- max(xy_resdev_post$x_upper, xy_resdev_post$y_upper)

      p <- ggplot2::ggplot(xy_resdev_post,
                           ggplot2::aes(y = .data$resdev_y,
                                        x = .data$resdev_x,
                                        xmin = .data$x_lower, xmax = .data$x_upper,
                                        ymin = .data$y_lower, ymax = .data$y_upper,
                                        label = .data$.label))
    } else {
      ulim <- max(xy_resdev_post$resdev_x, xy_resdev_post$resdev_y)

      p <- ggplot2::ggplot(xy_resdev_post,
                           ggplot2::aes(y = .data$resdev_y,
                                        x = .data$resdev_x,
                                        label = .data$.label))
    }

    if (dplyr::n_distinct(xy_resdev_post$Type) > 1) {
      p <- p +
        ggplot2::aes(colour = .data$Type) +
        ggplot2::scale_colour_viridis_d("")
    }

    p <- p +
      ggplot2::labs(x = "Residual Deviance (model 1)",
                    y = "Residual Deviance (model 2)") +
      ggplot2::coord_fixed(xlim = c(0, ulim), ylim = c(0, ulim)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60")

    if (show_uncertainty) {
      p <- p +
        # Have to layer up by data type by hand, due to the dual vertical and horizontal geoms
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "IPD"), orientation = "vertical")) +
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "IPD"), orientation = "horizontal")) +
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "AgD (arm-based)"), orientation = "vertical")) +
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "AgD (arm-based)"), orientation = "horizontal")) +
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "AgD (contrast-based)"), orientation = "vertical")) +
        do.call(tb_geom, args = purrr::list_modify(geom_dots, data = ~dplyr::filter(., .data$Type == "AgD (contrast-based)"), orientation = "horizontal"))
    } else {
      p <- p + ggplot2::geom_point(...)
    }

    p <- p + theme_multinma()

  } else {

    if (type == "resdev") { # Produce resdev plot

      if (show_uncertainty) {
        if (!rlang::is_string(stat))
          abort("`stat` should be a character string specifying the name of a ggdist stat")

        stat <- stringr::str_remove(stat, "^(stat_dist_|stat_|geom_)")

        tb_stat <- tryCatch(getExportedValue("ggdist", paste0("stat_", stat)),
                            error = function(err) {
                              abort(paste("`stat` should be a character string specifying the name of a ggdist stat:",
                                          err, sep = "\n"))
                            })

      }

      # Is a horizontal geom specified?
      horizontal <- orientation == "horizontal"

      if (horizontal) {
        resdev_post$.label <- forcats::fct_rev(resdev_post$.label)

        p <- ggplot2::ggplot(resdev_post,
                             ggplot2::aes(y = .data$.label,
                                          x = .data$resdev)) +
          ggplot2::labs(x = "Residual Deviance", y = "Data Point")

        if (dplyr::n_distinct(resdev_post$df) > 1) {
          if (dplyr::n_distinct(resdev_post$Type) > 1) {
            p <- p + ggplot2::facet_grid(Type~df_label, space = "free", scales = "free_y") +
              ggplot2::geom_vline(ggplot2::aes(xintercept = df), colour = "grey60")
          } else {
            p <- p + ggplot2::facet_grid(df_label~., space = "free", scales = "free_y") +
              ggplot2::geom_vline(ggplot2::aes(xintercept = df), colour = "grey60")
          }
        } else {
          p <- p + ggplot2::geom_vline(xintercept = 1, colour = "grey60")
          if (dplyr::n_distinct(resdev_post$Type) > 1) {
            p <- p + ggplot2::facet_grid(Type~., space = "free", scales = "free_y")
          }
        }

      } else {
        p <- ggplot2::ggplot(resdev_post,
                             ggplot2::aes(x = .data$.label,
                                          y = .data$resdev)) +
          ggplot2::labs(y = "Residual Deviance", x = "Data Point")

        if (dplyr::n_distinct(resdev_post$df) > 1) {
          if (dplyr::n_distinct(resdev_post$Type) > 1) {
            p <- p + ggplot2::facet_grid(df_label~Type, space = "free", scales = "free_x") +
              ggplot2::geom_hline(ggplot2::aes(yintercept = df), colour = "grey60")
          } else {
            p <- p + ggplot2::facet_grid(.~df_label, space = "free", scales = "free_x") +
              ggplot2::geom_hline(ggplot2::aes(yintercept = df), colour = "grey60")
          }
        } else {
          p <- p + ggplot2::geom_hline(yintercept = 1, colour = "grey60")
          if (dplyr::n_distinct(resdev_post$Type) > 1) {
            p <- p + ggplot2::facet_grid(.~Type, space = "free", scales = "free_x")
          }
        }
      }

      if (show_uncertainty) {
        p <- p + do.call(tb_stat,
                         args = rlang::dots_list(point_interval = ggdist::mean_qi,
                                                 ...,
                                                 orientation = orientation,
                                                 .homonyms = "last"))
      } else {
        p <- p + ggplot2::geom_point(...)
      }

      p <- p + theme_multinma()

      if (!horizontal) {
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
      } else {
        p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0, vjust = 0.5))
      }

    } else { # Leverage plot

      # Calculate signed sqrt resdev, add in leverages
      if (NROW(x$pointwise$ipd)) {
        resdev_post_ipd <- dplyr::right_join(resdev_post,
                                        dplyr::transmute(x$pointwise$ipd,
                                                         parameter = make_data_labels(.data$.study, .data$.trt),
                                                         .data$dic,
                                                         .data$leverage,
                                                         ssrd =
                                                           if (is.matrix(.data$observed)) {  # Ordered multinomial case
                                                            sign(rowSums(.data$observed[,-1] - .data$fitted[,-1], na.rm = TRUE)) * sqrt(.data$resdev)
                                                           } else {
                                                             sign(.data$observed - .data$fitted) * sqrt(.data$resdev)
                                                           }
                                        ),
                                        by = "parameter")
      }
      if (NROW(x$pointwise$agd_arm)) {
        resdev_post_aa <- dplyr::right_join(resdev_post,
                                        dplyr::transmute(x$pointwise$agd_arm,
                                                         parameter = make_data_labels(.data$.study, .data$.trt),
                                                         .data$dic,
                                                         .data$leverage,
                                                         ssrd =
                                                           if (is.matrix(.data$observed)) {  # Ordered multinomial case
                                                             sign(rowSums(.data$observed[,-1] - .data$fitted[,-1], na.rm = TRUE)) * sqrt(.data$resdev)
                                                           } else {
                                                             sign(.data$observed - .data$fitted) * sqrt(.data$resdev)
                                                           }
                                        ),
                                        by = "parameter")
      }
      if (NROW(x$pointwise$agd_contrast)) {
        resdev_post_ac <- dplyr::right_join(resdev_post,
                                        dplyr::transmute(x$pointwise$agd_contrast,
                                                         parameter = .data$.study,
                                                         .data$dic,
                                                         .data$leverage,
                                                         ssrd = purrr::map2_dbl(.data$observed,
                                                                                .data$fitted,
                                                                                ~sign(sum(.x - .y))) * sqrt(.data$resdev)
                                        ),
                                        by = "parameter")
      }

      resdev_post <- dplyr::bind_rows(if (NROW(x$pointwise$ipd)) resdev_post_ipd else tibble::tibble(),
                                      if (NROW(x$pointwise$agd_arm)) resdev_post_aa else tibble::tibble(),
                                      if (NROW(x$pointwise$agd_contrast)) resdev_post_ac else tibble::tibble())

      # Standardise sqrt residual deviances and leverages by df (e.g. for multi-arm trials)
      resdev_post$ssrd <- resdev_post$ssrd / sqrt(resdev_post$df)
      resdev_post$leverage <- resdev_post$leverage / resdev_post$df


      p <- ggplot2::ggplot(resdev_post,
                           ggplot2::aes(y = .data$leverage,
                                        x = .data$ssrd,
                                        label = .data$.label))

      if (dplyr::n_distinct(resdev_post$Type) > 1) {
        p <- p +
          ggplot2::aes(colour = .data$Type) +
          ggplot2::scale_colour_viridis_d("")
      }

      p <- p +
        ggplot2::labs(x = "Signed Square Root Residual Deviance",
                      y = "Leverage")

      rmax <- max(sqrt(resdev_post$resdev))

      if (!all(is.na(dic_contours))) {
        p <- p +
          purrr::map(dic_contours,
                     ~ggplot2::geom_function(fun = function(x, c) c - x^2, args = list(c = .),
                                            colour = "grey60", inherit.aes = FALSE,
                                            xlim = c(-rmax, rmax)*1.2)) +
          ggplot2::annotate("label", vjust = 0, label.size = 0,  colour = "grey60", fill = NA,
                            x = 0, y = dic_contours, label = paste0("DIC = ", dic_contours))
      }

      p <- p +
        ggplot2::geom_point(...) +
        ggplot2::coord_cartesian(
          ylim = c(min(0, resdev_post$leverage, na.rm = TRUE),
                   max(dic_contours, resdev_post$leverage, na.rm = TRUE)),
          xlim = c(-rmax, rmax)) +
        theme_multinma()

    }
  }

  return(p)
}
