#' Set up individual patient data
#'
#' @template args-data_common
#' @template args-data_arm
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @seealso [set_agd_arm()] for arm-based aggregate data, [set_agd_contrast()]
#'   for contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @examples
set_ipd <- function(data,
                    study,
                    trt,
                    y = NULL,
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

  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- dplyr::pull(data, {{ study }})
  if (any(is.na(.study))) abort("`study` cannot contain missing values")

  if (missing(trt)) abort("Specify `trt`")
  .trt <- dplyr::pull(data, {{ trt }})
  if (any(is.na(.trt))) abort("`trt` cannot contain missing values")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .r <- pull_non_null(data, enquo(r))
  .E <- pull_non_null(data, enquo(E))
  # .Surv <- ...

  o_continuous <- check_outcome_continuous(.y, with_se = FALSE)
  o_binary <- check_outcome_binary(.r, .E)
  # o_surv <- check_outcome_surv(.Surv)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = factor(.study),
    .trt = factor(.trt),
    .y = o_continuous$y,
    .r = o_binary$r,
    .E = o_binary$E#,
    # .Surv = o_surv$Surv,
  )

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = d,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study)),
    class = "nma_data")
  return(out)
}


#' Set up arm-based aggregate data
#'
#' @template args-data_common
#' @template args-data_arm
#' @param se column of `data` specifying the standard error for a continuous
#'   outcome
#' @param n column of `data` specifying Binomial outcome numerator
#'
#' @return An object of class [nma_data]
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

  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- dplyr::pull(data, {{ study }})
  if (any(is.na(.study))) abort("`study` cannot contain missing values")

  if (missing(trt)) abort("Specify `trt`")
  .trt <- dplyr::pull(data, {{ trt }})
  if (any(is.na(.trt))) abort("`trt` cannot contain missing values")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .se <- pull_non_null(data, enquo(se))
  .r <- pull_non_null(data, enquo(r))
  .n <- pull_non_null(data, enquo(n))
  .E <- pull_non_null(data, enquo(E))
  # .Surv <- ...

  o_continuous <- check_outcome_continuous(.y, .se, with_se = TRUE)
  o_count <- check_outcome_count(.r, .n, .E)
  # o_surv <- check_outcome_surv(.Surv)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = factor(.study),
    .trt = factor(.trt),
    .y = o_continuous$y,
    .se = o_continuous$se,
    .r = o_count$r,
    .n = o_count$n,
    .E = o_count$E#,
    # .Surv = o_surv$Surv,
    )

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = d,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study)),
    class = "nma_data")
  return(out)
}


#' Set up contrast-based aggregate data
#'
#' @template args-data_common
#' @param se column of `data` specifying the standard error for a continuous
#'   outcome
#' @param trt_b column of `data` specifying the reference/baseline treatment for
#'   each contrast
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_arm()] for
#'   arm-based aggregate data, and [combine_network()] for combining several
#'   data sources in one network.
#' @examples
set_agd_contrast <- function(data,
                             study,
                             trt, trt_b,
                             y = NULL, se = NULL) {

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


  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- dplyr::pull(data, {{ study }})
  if (any(is.na(.study))) abort("`study` cannot contain missing values")

  if (missing(trt)) abort("Specify `trt`")
  .trt <- dplyr::pull(data, {{ trt }})
  if (any(is.na(.trt))) abort("`trt` cannot contain missing values")

  if (missing(trt_b)) abort("Specify `trt_b`")
  .trt_b <- dplyr::pull(data, {{ trt_b }})
  if (any(is.na(.trt_b))) abort("`trt_b` cannot contain missing values")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .se <- pull_non_null(data, enquo(se))

  o_continuous <- check_outcome_continuous(.y, .se, with_se = TRUE)

  # Get all treatments
  trts <- sort(forcats::lvls_union(list(factor(.trt), factor(.trt_b))))

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = factor(.study),
    .trt = factor(.trt, levels = trts),
    .trt_b = factor(.trt_b, levels = trts),
    .y = o_continuous$y,
    .se = o_continuous$se)

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = d,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study)),
    class = "nma_data")
  return(out)
}


#' Combine multiple data sources into one network
#'
#' @param ... multiple data sources, as defined using the `set_*` functions
#' @param trt_ref reference treatment for the entire network, as a string (or
#'   coerced as such) referring to the levels of the treatment factor variable
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @seealso [set_ipd()], [set_agd_arm()], and [set_agd_contrast()] for defining
#'   different data sources
#' @examples
combine_network <- function(..., trt_ref, trt_refn) {
  s <- list(...)

  # Check that arguments all inherit from nma_data class
  if (!purrr::every(s, inherits, what = "nma_data")) {
    abort("Expecting to combine objects of class `nma_data`, created using set_* functions")
  }

  # Combine treatment code factor
  trts <- sort(forcats::lvls_union(purrr::map(s, "treatments")))
  if (!missing(trt_ref)) {
    if (! trt_ref %in% trts) {
      abort(sprintf("`trt_ref` does not match a treatment in the network.\nSuitable values are: %s",
                      ifelse(length(trts) <= 5,
                             paste0(trts, collapse = ", "),
                             paste0(paste0(trts[1:5], collapse = ", "), ", ..."))))
    }
    trts <- c(trt_ref, setdiff(trts, trt_ref))
  }

  # Check that no studies are duplicated between data sources
  all_studs <- purrr::flatten_chr(purrr::map(s, ~levels(.$studies)))
  if (anyDuplicated(all_studs)) {
    abort(sprintf("Studies with same label found in multiple data sources: %s",
                  paste0(unique(all_studs[duplicated(all_studs)]), collapse = ", ")))
  }

  # Combine study code factor
  studs <- sort(forcats::lvls_union(purrr::map(s, "studies")))

  # Get ipd
  ipd <- purrr::map_dfr(s, "ipd")
  if (!rlang::is_empty(ipd)) {
    ipd$.trt <- forcats::lvls_expand(ipd$.trt, trts)
    ipd$.study <- forcats::lvls_expand(ipd$.study, studs)
  }

  # Get agd_arm
  agd_arm <- purrr::map_dfr(s, "agd_arm")
  if (!rlang::is_empty(agd_arm)) {
    agd_arm$.trt <- forcats::lvls_expand(agd_arm$.trt, trts)
    agd_arm$.study <- forcats::lvls_expand(agd_arm$.study, studs)
  }

  # Get agd_contrast
  agd_contrast <- purrr::map_dfr(s, "agd_contrast")
  if (!rlang::is_empty(agd_contrast)) {
    agd_contrast$.trt <- forcats::lvls_expand(agd_contrast$.trt, trts)
    agd_contrast$.study <- forcats::lvls_expand(agd_contrast$.study, studs)
  }

  # Produce nma_data object
  out <- structure(
    list(agd_arm = agd_arm,
         agd_contrast = agd_contrast,
         ipd = ipd,
         treatments = factor(trts, levels = trts),
         studies = factor(studs)),
    class = "nma_data")
  return(out)
}

#' Pull non-null variables from data
#'
#' @param data data frame
#' @param var quosure (possibly NULL) for variable to pull
#'
#' @noRd
pull_non_null <- function(data, var) {
  var_null <- rlang::quo_is_missing(var) | rlang::quo_is_null(var)
  if (!var_null) return(dplyr::pull(data, {{ var }}))
  else return(NULL)
}


#' Check continuous outcomes
#'
#' @param y vector
#' @param se vector
#' @param with_se continuous outcome with SE?
#'
#' @noRd
check_outcome_continuous <- function(y, se = NULL, with_se = TRUE) {
  null_y <- is.null(y)
  null_se <- is.null(se)

  if (with_se) {
    if (!null_y & !null_se) {
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
  } else {
    if (!null_y) {
      if (!is.numeric(y)) abort("Continuous outcome `y` must be numeric")
    } else {
      y <- NA_real_
    }
    return(list(y = y))
  }
}

#' Check count outcomes
#'
#' @param r vector
#' @param n vector
#' @param E vector
#'
#' @noRd
check_outcome_count <- function(r, n, E) {
  null_r <- is.null(r)
  null_n <- is.null(n)
  null_E <- is.null(E)

  if (!null_r & !null_n) {
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
        if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
        if (any(E <= 0)) abort("Time at risk `E` must be positive")
      }
  } else {
    E <- NA_real_
  }
  return(list(r = r, n = n, E = E))
}

#' Check binary outcomes
#'
#' @param r vector
#' @param E vector
#'
#' @noRd
check_outcome_binary <- function(r, E) {
  null_r <- is.null(r)
  null_E <- is.null(E)

  if (!null_r) {
    if (!is.numeric(r)) abort("Binary outcome `r` must be numeric")
    if (any(! r %in% c(0, 1))) abort("Binary outcome `r` must equal 0 or 1")
  } else {
    r <- NA_integer_
  }
  if (!null_E) {
    if (null_r) {
      warn("Ignoring `E` without `r`")
      E <- NA_real_
    } else {
      if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
      if (any(E <= 0)) abort("Time at risk `E` must be positive")
    }
  } else {
    E <- NA_real_
  }
  return(list(r = r, E = E))
}
