#' Set up individual patient data
#'
#' @templateVar args-data data, trt, study, y, se, r, E, surv
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
                    surv) {
  # ...

  # Produce nma_data_ipd object
  out <- structure(list(), class = c("nma_data_ipd", "nma_data"))
  return(out)
}


#' Set up arm-based aggregate data
#'
#' @templateVar args-data data, trt, study, y, se, r, n, E, surv
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
                        surv) {
  # ...

  # Produce nma_data_agd_arm object
  out <- structure(list(), class = c("nma_data_agd_arm", "nma_data"))
  return(out)
}


#' Set up contrast-based aggregate data
#'
#' @templateVar args-data data, trt, study, y, se
#' @param trt_b
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
#' @param ...
#' @param trt_ref
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
