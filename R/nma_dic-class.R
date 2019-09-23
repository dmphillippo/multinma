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
#'   \item{`dbar`}{The total residual deviance}
#'   \item{`pointwise`}{A list of data frames containing the pointwise
#'   contributions for the IPD and AgD.}
#'   }
#'
NULL
