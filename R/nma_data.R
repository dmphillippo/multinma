#' Set up individual patient data
#'
#' @template args-data_common
#' @template args-data_arm
#'
#' @return
#' @export
#'
#' @seealso [set_agd_arm()] for arm-based aggregate data, [set_agd_contrast()]
#'   for contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @examples
set_ipd <- function(data,
                    study,
                    trt,
                    y = NULL, se = NULL,
                    r = NULL, E = NULL,
                    Surv = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  # ...

  dat <- tibble::tibble(
    .study = factor(study),
    .trt = factor(trt)
  )

  # Produce nma_data_ipd object
  out <- structure(list(), class = "nma_data")
  return(out)
}


#' Set up arm-based aggregate data
#'
#' @template args-data_common
#' @template args-data_arm
#' @param n column of `data` specifying Binomial outcome numerator
#'
#' @return
#' @export
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_contrast()] for
#'   contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @examples
set_agd_arm <- function(data,
                        study,
                        trt,
                        y = NULL, se = NULL,
                        r = NULL, n = NULL, E = NULL,
                        Surv = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  # Pull single columns
  .study <- dplyr::pull(data, {{ study }})
  .trt <- dplyr::pull(data, {{ trt }})

  # Allow any number of outcomes to be given, but check that all components are present
  o_continuous <- check_outcome_continuous(data, enquo(y), enquo(se))
  o_count <- check_outcome_count(data, enquo(r), enquo(n), enquo(E))
  # o_surv <- check_outcome_surv(data, enquo(Surv))

  # Create tibble in standard format
  dat <- tibble::tibble(
    .study = factor(.study),
    .trt = factor(.trt),
    .y = o_continuous$y,
    .se = o_continuous$se,
    .r = o_count$r,
    .n = o_count$n,
    .E = o_count$E#,
    # .Surv = o_surv$Surv
  )

  # Produce nma_data_agd_arm object
  out <- structure(
    list(agd_arm = dat,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = NA,
         studies = NA),
    class = "nma_data")
  return(out)
}


#' Set up contrast-based aggregate data
#'
#' @template args-data_common
#' @param trt_b column of `data` specifying the reference/baseline treatment for
#'   each contrast
#'
#' @return
#' @export
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_arm()] for
#'   arm-based aggregate data, and [combine_network()] for combining several
#'   data sources in one network.
#' @examples
set_agd_contrast <- function(data,
                             study,
                             trt, trt_b,
                             y, se) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  # ...

  dat <- tibble::tibble(
    .study = factor(study),
    .trt = factor(trt)
  )

  # Produce nma_data_agd_contrast object
  out <- structure(list(), class = "nma_data")
  return(out)
}


#' Combine multiple data sources into one network
#'
#' @param ... multiple data sources, as defined using the `set_*` functions
#' @param trt_ref reference treatment for the entire network (string, factor, or
#'   integer)
#'
#' @return
#' @export
#'
#' @seealso [set_ipd()], [set_agd_arm()], and [set_agd_contrast()] for defining
#'   different data sources
#' @examples
combine_network <- function(..., trt_ref) {
  # Check that arguments all inherit from nma_data class

  # Check that all data sources use same reference treatment, warn otherwise to
  # use trt_ref

  # Combine treatment code factor

  # Combine study code factor

  # Produce nma_data object
  out <- structure(list(), class = "nma_data")
  return(out)
}

#' Check continuous outcomes
#'
#' @param data data frame
#' @param y quosure
#' @param se quosure
#'
#' @noRd
check_outcome_continuous <- function(data, y, se) {
  # Check NULL instead of missing(), since "missingness" is not passed from
  # calling function
  null_y <- rlang::quo_is_null(y)
  null_se <- rlang::quo_is_null(se)

  if (!null_y & !null_se) {
    y <- dplyr::pull(data, {{ y }})
    se <- dplyr::pull(data, {{ se }})

    if (!is.numeric(y)) abort("Continuous outcome `y` must be numeric")
    if (!is.numeric(se)) abort("Standard error `se` must be numeric")
    if (any(se <= 0)) abort("Standard errors must be positive")
  } else {
    if (!null_y) abort("Specify standard error `se` for continuous outcome `y`")
    if (!null_se) warn("Ignoring standard error `se` without continuous outcome `y`")
    y <- NA_real_
    se <- NA_real_
  }
  return(list(y = y, se = se))
}

#' Check count outcomes
#'
#' @param data data frame
#' @param r quosure
#' @param n quosure
#' @param E quosure
#'
#' @noRd
check_outcome_count <- function(data, r, n, E) {
  # Check NULL instead of missing(), since "missingness" is not passed from
  # calling function
  null_r <- rlang::quo_is_null(r)
  null_n <- rlang::quo_is_null(n)
  null_E <- rlang::quo_is_null(E)

  if (!null_r & !null_n) {
    r <- dplyr::pull(data, {{ r }})
    n <- dplyr::pull(data, {{ n }})

    if (!is.numeric(r)) abort("Count outcome `r` must be numeric")
    if (!is.numeric(n)) abort("Denominator `n` must be numeric")
    if (any(r != trunc(r))) abort("Count outcome `r` must be integer-valued")
    if (any(n != trunc(n))) abort("Denominator `n` must be integer-valued")
    if (any(n <= 0)) abort("Denominator `n` must be greater than zero")
    if (any(n < r | r < 0)) abort("Count outcome `r` must be between 0 and `n`")
  } else if (!null_n & null_r) {
    warn("Ignoring `n` without `r`")

    r <- NA_integer_
    n <- NA_integer_
  } else if (null_n & !null_r) {
    abort("Specify denominator `n` for count outcome `r`")
  } else {
    r <- NA_integer_
    n <- NA_integer_
  }
  if (!null_E) {
      if (null_n | null_r) {
        warn("Ignoring `E` without `r` or `n`")
        E <- NA_real_
      } else {
        E <- dplyr::pull(data, {{ E }})
        if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
        if (any(E <= 0)) abort("Time at risk `E` must be positive")
      }
  } else {
    E <- NA_real_
  }
  return(list(r = r, n = n, E = E))
}
