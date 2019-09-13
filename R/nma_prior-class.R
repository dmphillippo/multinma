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
#'   \item{`...`}{Parameters of the distribution}
#'   }
#'
#' The distribution parameters, specified as named components in `...`, match
#' those in the constructor functions (see [priors]).
NULL

#' @export
#' @noRd
print.nma_prior <- function(x, ...) {
  p <- purrr::list_modify(x, dist = purrr::zap())
  p <- p[!is.na(p)]
  cglue("A {x$dist} prior distribution: {paste(names(p), p, sep = ' = ', collapse = ', ')}.")
  invisible(x)
}
