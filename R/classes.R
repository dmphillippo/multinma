#' The nma_data class
#'
#' The `nma_data` class contains the data for a NMA in a standard format,
#' created using the functions [set_ipd()], [set_agd_arm()],
#' [set_agd_contrast()], or [combine()].
#'
#' @rdname nma_data-class
#' @name nma_data-class
#'
#' @details Objects of class `nma_data` have the following components:
#'   \describe{
#'   \item{`agd_arm`}{data from studies with aggregate data (arm format)}
#'   \item{`agd_contrast`}{data from studies with aggregate data (contrast
#'   format)}
#'   \item{`ipd`}{data from studies with individual patient data}
#'   \item{`treatments`}{treatment coding factor for entire network}
#'   \item{`studies`}{study coding factor for entire network}
#'   }
#'
#' The `agd_arm`, `agd_contrast`, and `ipd` components are
#' tibbles with the following columns:
#'   \describe{
#'   \item{`.study`}{study (as factor)}
#'   \item{`.trt`}{treatment (as factor)}
#'   \item{`.trt_b`}{baseline treatment for contrast data (as factor)}
#'   \item{`.y`}{outcome (continuous)}
#'   \item{`.se`}{standard error (continuous)}
#'   \item{`.r`}{event count (discrete)}
#'   \item{`.n`}{total number of individuals (discrete, `agd_arm` only)}
#'   \item{`.E`}{time at risk (discrete)}
#'   \item{`.surv`}{event/censoring time, of type `Surv` (time-to-event)}
#'   \item{`...`}{other columns (typically covariates) from the original data
#'   frame}
#'   }
#'
#' The classes `nma_data_ipd`, `nma_data_agd_arm`, and `nma_data_agd_contrast`
#' inherit from the `nma_data` class, with only the corresponding data
#' components filled.
#'
NULL
