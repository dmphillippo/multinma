#' Set up individual patient data
#'
#' @template args-data_common
#' @template args-data_rE
#' @template args-data_Surv
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

  check_outcome_continuous(.y, with_se = FALSE)
  check_outcome_binary(.r, .E)
  # check_outcome_surv(.Surv)

  o_type <- get_outcome_type(y = .y, se = NULL,
                             r = .r, n = NULL, E = .E)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt)
  )

  if (o_type == "continuous") {
    d <- tibble::add_column(d, .y = .y)
  } else if (o_type == "binary") {
    d <- tibble::add_column(d, .r = .r)
  } else if (o_type == "rate") {
    d <- tibble::add_column(d, .r = .r, .E = .E)
  }

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = d,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = NA, agd_contrast = NA, ipd = o_type)),
    class = "nma_data")
  return(out)
}


#' Set up arm-based aggregate data
#'
#' @template args-data_common
#' @template args-data_se
#' @template args-data_rE
#' @template args-data_Surv
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

  check_outcome_continuous(.y, .se, with_se = TRUE)
  check_outcome_count(.r, .n, .E)
  # check_outcome_surv(.Surv)

  o_type <- get_outcome_type(y = .y, se = .se,
                             r = .r, n = .n, E = .E)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt)
  )

  if (o_type == "continuous") {
    d <- tibble::add_column(d, .y = .y, .se = .se)
  } else if (o_type == "count") {
    d <- tibble::add_column(d, .r = .r, .n = .n)
  } else if (o_type == "rate") {
    d <- tibble::add_column(d, .r = .r, .E = .E)
  }

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = d,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = o_type, agd_contrast = NA, ipd = NA)),
    class = "nma_data")
  return(out)
}


#' Set up contrast-based aggregate data
#'
#' @template args-data_common
#' @template args-data_se
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

  check_outcome_continuous(.y, .se, with_se = TRUE)

  o_type <- get_outcome_type(y = .y, se = .se,
                             r = NULL, n = NULL, E = NULL)

  # Get all treatments
  trts <- stringr::str_sort(forcats::lvls_union(list(nfactor(.trt), nfactor(.trt_b))), numeric = TRUE)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = factor(.trt, levels = trts),
    .trt_b = factor(.trt_b, levels = trts),
    .y = .y,
    .se = .se)

  d <- dplyr::bind_cols(d, data)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = d,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = NA, agd_contrast = o_type, ipd = NA)),
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
combine_network <- function(..., trt_ref) {
  s <- list(...)

  # Check that arguments all inherit from nma_data class
  if (!purrr::every(s, inherits, what = "nma_data")) {
    abort("Expecting to combine objects of class `nma_data`, created using set_* functions")
  }

  # Combine treatment code factor
  trts <- stringr::str_sort(forcats::lvls_union(purrr::map(s, "treatments")), numeric = TRUE)
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
  studs <- stringr::str_sort(forcats::lvls_union(purrr::map(s, "studies")), numeric = TRUE)

  # Get ipd
  ipd <- purrr::map(s, "ipd")
  if (!rlang::is_empty(ipd)) {
    for (j in 1:length(ipd)) {
      if (rlang::is_empty(ipd[[j]])) next
      ipd[[j]]$.trt <- forcats::lvls_expand(ipd[[j]]$.trt, trts)
      ipd[[j]]$.study <- forcats::lvls_expand(ipd[[j]]$.study, studs)
    }
  }
  ipd <- dplyr::bind_rows(ipd)

  # Get agd_arm
  agd_arm <- purrr::map(s, "agd_arm")
  if (!rlang::is_empty(agd_arm)) {
    for (j in 1:length(agd_arm)) {
      if (rlang::is_empty(agd_arm[[j]])) next
      agd_arm[[j]]$.trt <- forcats::lvls_expand(agd_arm[[j]]$.trt, trts)
      agd_arm[[j]]$.study <- forcats::lvls_expand(agd_arm[[j]]$.study, studs)
    }
  }
  agd_arm <- dplyr::bind_rows(agd_arm)

  # Get agd_contrast
  agd_contrast <- purrr::map(s, "agd_contrast")
  if (!rlang::is_empty(agd_contrast)) {
    for (j in 1:length(agd_contrast)) {
      if (rlang::is_empty(agd_contrast[[j]])) next
      agd_contrast[[j]]$.trt <- forcats::lvls_expand(agd_contrast[[j]]$.trt, trts)
      agd_contrast[[j]]$.study <- forcats::lvls_expand(agd_contrast[[j]]$.study, studs)
    }
  }
  agd_contrast <- dplyr::bind_rows(agd_contrast)

  # Get outcome type
  o_ipd <- unique(purrr::map_chr(purrr::map(s, "outcome"), "ipd"))
  o_ipd <- o_ipd[!is.na(o_ipd)]
  if (length(o_ipd) > 1) abort("Multiple outcome types present in IPD.")
  if (length(o_ipd) == 0) o_ipd <- NA

  o_agd_arm <- unique(purrr::map_chr(purrr::map(s, "outcome"), "agd_arm"))
  o_agd_arm <- o_agd_arm[!is.na(o_agd_arm)]
  if (length(o_agd_arm) > 1) abort("Multiple outcome types present in AgD (arm-based).")
  if (length(o_agd_arm) == 0) o_agd_arm <- NA

  o_agd_contrast <- unique(purrr::map_chr(purrr::map(s, "outcome"), "agd_contrast"))
  o_agd_contrast <- o_agd_contrast[!is.na(o_agd_contrast)]
  if (length(o_agd_contrast) > 1) abort("Multiple outcome types present in AgD (contrast-based).")
  if (length(o_agd_contrast) == 0) o_agd_contrast <- NA

  outcome <- list(agd_arm = o_agd_arm,
                  agd_contrast = o_agd_contrast,
                  ipd = o_ipd)

  # Check outcome combination
  check_outcome_combination(outcome)

  # Produce nma_data object
  out <- structure(
    list(agd_arm = agd_arm,
         agd_contrast = agd_contrast,
         ipd = ipd,
         treatments = factor(trts, levels = trts),
         studies = factor(studs, levels = studs),
         outcome = outcome),
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

#' Get outcome type
#'
#' Determines outcome type based on which inputs are NA
#'
#' @noRd
get_outcome_type <- function(y, se, r, n, E) {
  o <- c()
  if (!is.null(y)) o <- c(o, "continuous")
  if (!is.null(r)) {
    if (!is.null(E)) o <- c(o, "rate")
    if (!is.null(n)) o <- c(o, "count")
    if (is.null(n) && is.null(E)) o <- c(o, "binary")
  }
  if (length(o) == 0) abort("Please specify one and only one outcome.")
  if (length(o) > 1) abort(glue::glue("Please specify one and only one outcome, instead of ",
                                      glue::glue_collapse(o, sep = ", ", last = " and "), "."))

  return(o)
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
    if (!null_y && !null_se) {
      if (!is.numeric(y)) abort("Continuous outcome `y` must be numeric")
      if (!is.numeric(se)) abort("Standard error `se` must be numeric")
      if (any(is.na(y))) abort("Continuous outcome `y` contains missing values")
      if (any(is.na(se))) abort("Standard error `se` contains missing values")
      if (any(se <= 0)) abort("Standard errors must be positive")
    } else {
      if (!null_y) abort("Specify standard error `se` for continuous outcome `y`")
      if (!null_se) abort("Specify continuous outcome `y`")
    }
    invisible(list(y = y, se = se))
  } else {
    if (!null_y) {
      if (any(is.na(y))) abort("Continuous outcome `y` contains missing values")
      if (!is.numeric(y)) abort("Continuous outcome `y` must be numeric")
    }
    invisible(list(y = y))
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

  if (!null_n) {
    if (!is.numeric(n)) abort("Denominator `n` must be numeric")
    if (any(is.na(n))) abort("Denominator `n` contains missing values")
    if (any(n != trunc(n))) abort("Denominator `n` must be integer-valued")
    if (any(n <= 0)) abort("Denominator `n` must be greater than zero")
    if (null_r) abort("Specify outcome count `r`.")
  }

  if (!null_E) {
    if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
    if (any(is.na(E))) abort("Time at risk `E` contains missing values")
    if (any(E <= 0)) abort("Time at risk `E` must be positive")
    if (null_r) abort("Specify outcome count `r`.")
  }

  if (!null_r) {
    if (null_n && null_E) abort("Specify denominator `n` (count outcome) or time at risk `E` (rate outcome)")
    if (!is.numeric(r)) abort("Outcome count `r` must be numeric")
    if (any(is.na(r))) abort("Outcome count `r` contains missing values")
    if (any(r != trunc(r))) abort("Outcome count `r` must be integer-valued")
    if (!null_n && any(n < r | r < 0)) abort("Count outcome `r` must be between 0 and `n`")
    if (!null_E && any(r < 0)) abort("Rate outcome count `r` must be non-negative")
  }

  invisible(list(r = r, n = n, E = E))
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

  if (!null_E) {
    if (null_r) {
      abort("Specify count `r` for rate outcome")
    } else {
      if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
      if (any(is.na(E))) abort("Time at risk `E` contains missing values")
      if (any(E <= 0)) abort("Time at risk `E` must be positive")
      if (!is.numeric(r)) abort("Rate outcome count `r` must be numeric")
      if (any(is.na(r))) abort("Rate outcome count `r` contains missing values")
      if (any(r != trunc(r))) abort("Rate outcome count `r` must be non-negative integer")
      if (any(r < 0)) abort("Rate outcome count `r` must be non-negative integer")
    }
  } else if (!null_r) {
    if (!is.numeric(r)) abort("Binary outcome `r` must be numeric")
    if (any(is.na(r))) abort("Binary outcome `r` contains missing values")
    if (any(! r %in% c(0, 1))) abort("Binary outcome `r` must equal 0 or 1")
  }

  invisible(list(r = r, E = E))
}

#' Check valid outcome combination across data sources
#'
#' @param outcomes outcome list, see nma_data-class
#'
#' @noRd
check_outcome_combination <- function(outcomes) {
  valid <- list(
    list(agd_arm = c("count", NA),
         agd_contrast = c("continuous", NA),
         ipd = c("binary", NA)),
    list(agd_arm = c("rate", NA),
         agd_contrast = c("continuous", NA),
         ipd = c("rate", NA)),
    list(agd_arm = c("continuous", NA),
         agd_contrast = c("continuous", NA),
         ipd = c("continuous", NA))
  )

  if (!any(purrr::map_lgl(valid,
                 ~all(c(outcomes$agd_arm %in% .$agd_arm,
                        outcomes$agd_contrast %in% .$agd_contrast,
                        outcomes$ipd %in% .$ipd))))) {
    rlang::abort(glue::glue("Combining ",
                     glue::glue_collapse(outcomes[!is.na(outcomes)], sep = ', ', last = ' and '),
                     " outcomes is not supported."))
  }
}

#' Check for IPD and AgD in network
#'
#' @param network nma_data object
#'
#' @return logical TRUE/FALSE
#' @noRd
has_ipd <- function(network) {
  if (!inherits(network, "nma_data")) abort("Not nma_data object.")
  return(!rlang::is_empty(network$ipd))
}

has_agd_arm <- function(network) {
  if (!inherits(network, "nma_data")) abort("Not nma_data object.")
  return(!rlang::is_empty(network$agd_arm))
}

has_agd_contrast <- function(network) {
  if (!inherits(network, "nma_data")) abort("Not nma_data object.")
  return(!rlang::is_empty(network$agd_contrast))
}

#' Natural-order factors
#'
#' Produces factors with levels in natural sort order (i.e. 1 5 10 not 1 10 5)
#'
#' @noRd
nfactor <- function(x, ..., numeric = TRUE) {
  return(factor(x, levels = stringr::str_sort(unique(x), numeric = numeric), ...))
}
