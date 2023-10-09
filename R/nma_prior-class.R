#' The nma_prior class
#'
#' The `nma_prior` class is used to specify prior distributions.
#'
#' @rdname nma_prior-class
#' @name nma_prior-class
#' @aliases nma_prior
#'
#' @details Objects of class `nma_prior` have the following components:
#'   \describe{
#'   \item{`dist`}{Distribution name}
#'   \item{`fun`}{Name of constructor function, as string (e.g. `"normal"`)}
#'   \item{`...`}{Parameters of the distribution}
#'   }
#'
#' The distribution parameters, specified as named components in `...`, match
#' those in the constructor functions (see [priors]).
NULL

#' @export
#' @noRd
print.nma_prior <- function(x, ...) {
  if (x$dist == "flat (implicit)") {
    cglue("An implicit flat prior distribution over the entire parameter support.")
  } else {
    p <- purrr::list_modify(x, dist = purrr::zap(), fun = purrr::zap())
    p <- p[!is.na(p)]
    cglue("A{if (stringr::str_starts(x$dist, '[aeiouAEIOU]')) 'n' else ''} {x$dist} prior distribution: {paste(names(p), p, sep = ' = ', collapse = ', ')}.")
  }
  invisible(x)
}


#' Summary of prior distributions
#'
#' Print a summary of prior distribution details.
#'
#' @param object Prior distribution as a `nma_prior` object
#' @param ... Additional arguments, not used
#' @param probs Numeric vector of probabilities to calculate prior intervals
#' @param digits Number of digits to display
#' @param trunc Optional numeric vector of length 2, giving the truncation
#'   limits of the prior distribution. Useful if a real-valued prior is assigned
#'   to a positive-valued parameter, then `trunc = c(0, Inf)` will give the
#'   correct prior intervals. By default, truncation is not used.
#'
#' @return A data frame is returned invisibly, giving the prior intervals
#' @export
#'
#' @examples
#' summary(normal(location = 0, scale = 1))
#' summary(half_normal(scale = 1))
#' summary(log_normal(location = -3.93, scale = 1.51))
#'
#' # Truncation limits may be set, for example to restrict a prior to positive values
#' summary(normal(location = 0.5, scale = 1), trunc = c(0, Inf))
#'
summary.nma_prior <- function(object, ..., probs = c(0.5, 0.95), digits = 2, trunc = NULL) {
  check_probs(probs)
  if (!rlang::is_scalar_integerish(digits, finite = TRUE))
    abort("`digits` must be a single integer.")
  if (!is.null(trunc) && !rlang::is_double(trunc, n = 2) || any(is.na(trunc)))
    abort("`trunc` must be a length 2 numeric vector of truncation limits.")

  if (object$dist == "flat (implicit)") {
    if (is.null(trunc)) trunc <- c(-Inf, Inf)
    cglue("A flat prior on the parameter support {round(trunc[1], digits)} to {round(trunc[2], digits)}.")
    invisible(tibble::tibble(probs = 1, lower = trunc[1], upper = trunc[2]))
  } else {
    prior <- get_tidy_prior(object, trunc = trunc, ...) %>%
      tidyr::expand_grid(probs = probs) %>%
      dplyr::group_by(.data[["dist"]], .data[["probs"]]) %>%
      {if (stringr::str_starts(object$dist, "half-|Exponential")) {
        dplyr::summarise(., qfun = paste0("q", .data[["dist"]]),
                         lower = 0,
                         upper = do.call(.data[["qfun"]], args = rlang::list2(p = .data[["probs"]], !!! .[["args"]][[1]])))
      } else {
        dplyr::summarise(., qfun = paste0("q", .data[["dist"]]),
                         lower = do.call(.data[["qfun"]], args = rlang::list2(p = (1 - .data[["probs"]]) / 2, !!! .[["args"]][[1]])),
                         upper = do.call(.data[["qfun"]], args = rlang::list2(p = 1 - (1 - .data[["probs"]]) / 2, !!! .[["args"]][[1]])))
      }}

    print(object)
    cglue("{prior$probs*100}% of the prior density lies between {round(prior$lower, digits)} and {round(prior$upper, digits)}.")

    invisible(prior[ , c("probs", "lower", "upper")])
  }
}
