#' Set up individual patient data
#'
#' @template args-data_common
#' @template args-data_arm
#'
#' @return
#' @export
#'
#' @examples
set_ipd <- function(data,
                    trt,
                    study,
                    y, se,
                    r, E,
                    Surv) {
  # ...

  # Produce nma_data_ipd object
  out <- structure(list(), class = c("nma_data_ipd", "nma_data"))
  return(out)
}


#' Set up arm-based aggregate data
#'
#' @template args-data_common
#' @template args-data_arm
#' @param n column of \code{data} specifying Binomial outcome numerator
#'
#' @return
#' @export
#'
#' @examples
set_agd_arm <- function(data,
                        trt,
                        study,
                        y, se,
                        r, n, E,
                        Surv) {
  # ...

  # Produce nma_data_agd_arm object
  out <- structure(list(), class = c("nma_data_agd_arm", "nma_data"))
  return(out)
}


#' Set up contrast-based aggregate data
#'
#' @template args-data_common
#' @param trt_b column of \code{data} specifying the reference/baseline
#'   treatment for each contrast
#'
#' @return
#' @export
#'
#' @examples
set_agd_contrast <- function(data,
                             trt, trt_b,
                             study,
                             y, se) {
  # ...

  # Produce nma_data_agd_contrast object
  out <- structure(list(), class = c("nma_data_agd_contrast", "nma_data"))
  return(out)
}


#' Combine multiple data sources into one network
#'
#' @param ... multiple data sources, as defined using the \code{set_*} functions
#' @param trt_ref reference treatment for the entire network (string, factor, or
#'   integer)
#'
#' @return
#' @export
#'
#' @examples
combine <- function(..., trt_ref) {
  # Check that arguments all inherit from nma_data class

  # Check that all data sources use same reference treatment, warn otherwise to
  # use trt_ref

  # Combine treatment code factor

  # Combine study code factor

  # Produce nma_data object
  out <- structure(list(), class = "nma_data")
  return(out)
}
