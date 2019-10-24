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
