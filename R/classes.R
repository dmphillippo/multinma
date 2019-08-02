#' The nma_data class
#'
#' The \code{nma_data} class contains the data for a NMA in a standard format,
#' created using the functions \code{set_ipd}, \code{set_agd_arm},
#' \code{set_agd_contrast}, or \code{combine}.
#'
#' @rdname nma_data-class
#' @name nma_data-class
#'
#' @details Objects of class \code{nma_data} have the following components:
#'   \describe{
#'   \item{\code{agd_arm}}{data from studies with aggregate data (arm format)}
#'   \item{\code{agd_contrast}}{data from studies with aggregate data (contrast
#'   format)}
#'   \item{\code{ipd}}{data from studies with individual patient data}
#'   \item{\code{treatments}}{treatment coding factor for entire network}
#'   \item{\code{studies}}{study coding factor for entire network}
#'   }
#'
#' The \code{agd_arm}, \code{agd_contrast}, and \code{ipd} components are
#' tibbles with the following columns:
#'   \describe{
#'   \item{\code{.study}}{study (as factor)}
#'   \item{\code{.trt}}{treatment (as factor)}
#'   \item{\code{.trt_b}}{baseline treatment for contrast data (as factor)}
#'   \item{\code{.y}}{outcome (continuous)}
#'   \item{\code{.se}}{standard error (continuous)}
#'   \item{\code{.r}}{event count (discrete)}
#'   \item{\code{.n}}{total number of individuals (discrete, \code{agd_arm}
#'   only)}
#'   \item{\code{.E}}{time at risk (discrete)}
#'   \item{\code{.surv}}{event/censoring time, of type \code{Surv}
#'   (time-to-event)}
#'   \item{\code{...}}{other columns (typically covariates) from the original
#'   data frame}
#'   }
#'
#' The classes \code{nma_data_ipd}, \code{nma_data_agd_arm}, and
#' \code{nma_data_agd_contrast} inherit from the \code{nma_data} class, with
#' only the corresponding data components filled.
#'
NULL
