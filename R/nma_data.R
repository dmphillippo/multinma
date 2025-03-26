#' Set up individual patient data
#'
#' Set up a network containing individual patient data (IPD). Multiple data
#' sources may be combined once created using [combine_network()].
#'
#' @template args-data_common
#' @template args-data_y
# #' @template args-data_rE
#' @param r column of `data` specifying a binary outcome or Poisson outcome count
#' @param E column of `data` specifying the total time at risk for Poisson
#'   outcomes
#' @template args-data_Surv
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @template args-details_trt_ref
#' @template args-details_mutate
#'
#' @seealso [set_agd_arm()] for arm-based aggregate data, [set_agd_contrast()]
#'   for contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @template seealso_nma_data
#' @examples
#' # Set up network of plaque psoriasis IPD
#' head(plaque_psoriasis_ipd)
#'
#' pso_net <- set_ipd(plaque_psoriasis_ipd,
#'                    study = studyc,
#'                    trt = trtc,
#'                    r = pasi75)
#'
#' # Print network details
#' pso_net
#'
#' # Plot network
#' plot(pso_net)
#'
#' # Setting a different reference treatment
#' set_ipd(plaque_psoriasis_ipd,
#'         study = studyc,
#'         trt = trtc,
#'         r = pasi75,
#'         trt_ref = "PBO")

set_ipd <- function(data,
                    study,
                    trt,
                    y = NULL,
                    r = NULL, E = NULL,
                    Surv = NULL,
                    trt_ref = NULL,
                    trt_class = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             classes = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- pull_non_null(data, enquo(study))
  if (is.null(.study)) abort("`study` cannot be NULL")
  check_study(.study)

  if (is.factor(.study)) {
    study_original_levels <- levels(.study)
    .study <- forcats::fct_drop(.study)
  } else {
    study_original_levels <- NULL
  }

  if (missing(trt)) abort("Specify `trt`")
  .trt <- pull_non_null(data, enquo(trt))
  if (is.null(.trt)) abort("`trt` cannot be NULL")
  check_trt(.trt)

  if (is.factor(.trt)) {
    trt_original_levels <- levels(.trt)
    .trt <- forcats::fct_drop(.trt)
  } else {
    trt_original_levels <- NULL
  }

  # Treatment classes
  .trtclass <- pull_non_null(data, enquo(trt_class))
  if (!is.null(.trtclass)) {
    check_trt_class(.trtclass, .trt)

    if (is.factor(.trtclass)) {
      trtclass_original_levels <- levels(.trtclass)
      .trtclass <- forcats::fct_drop(.trtclass)
    } else {
      trtclass_original_levels <- NULL
    }
  }

  if (!is.null(trt_ref) && length(trt_ref) > 1) abort("`trt_ref` must be length 1.")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .r <- pull_non_null(data, enquo(r))
  .E <- pull_non_null(data, enquo(E))
  .Surv <- pull_non_null(data, enquo(Surv))

  check_outcome_continuous(.y, with_se = FALSE)

  if (!is.null(.r) && inherits(.r, c("multi_ordered", "multi_competing"))) {
    if (inherits(.r, "multi_competing")) abort("Competing multinomial outcomes are not yet supported.")

    # Most checks are carried out by multi(), some additional checks here specific to IPD
    if (any(! .r[!is.na(.r)] %in% c(0, 1))) abort("Multinomial outcome `r` must equal 0 or 1")
    if (any(r_zero_rows <- rowSums(.r, na.rm = TRUE) == 0))
      abort(glue::glue("Individual{if (sum(r_zero_rows) > 1) 's' else ''} without outcomes in any category, ",
                       "row{if (sum(r_zero_rows) > 1) 's' else ''} ",
                       glue::glue_collapse(which(r_zero_rows), sep = ", ", last = " and "), "."))
    if (any(r_multi_rows <- rowSums(.r, na.rm = TRUE) > 1))
      abort(glue::glue("Individual{if (sum(r_multi_rows) > 1) 's' else ''} with outcomes in more than one category, ",
                       "row{if (sum(r_multi_rows) > 1) 's' else ''} ",
                       glue::glue_collapse(which(r_multi_rows), sep = ", ", last = " and "), "."))
  } else {
    check_outcome_binary(.r, .E)
  }

  check_outcome_survival(.Surv)

  o_type <- get_outcome_type(y = .y, se = NULL,
                             r = .r, n = NULL, E = .E,
                             Surv = .Surv)

  # Check for single-arm studies
  single_arm_studies <- tibble::tibble(.study, .trt) %>%
    dplyr::distinct(.data$.study, .data$.trt) %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::pull(.data$.study)

  if (length(single_arm_studies)) {
    if (o_type == "survival") {
      inform(glue::glue("Single-arm stud{if (length(single_arm_studies) > 1) 'ies' else 'y'} present in the network: ",
                        glue::glue_collapse(glue::double_quote(as.character(single_arm_studies)), sep = ", ", last = " and "), "."))
    } else {
      abort(glue::glue("Single-arm studies are not supported: issue with stud{if (length(single_arm_studies) > 1) 'ies' else 'y'} ",
                       glue::glue_collapse(glue::double_quote(single_arm_studies), sep = ", ", last = " and "), "."))
    }
  }

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt)
  )

  if (!is.null(trt_ref)) {
    trt_ref <- as.character(trt_ref)
    lvls_trt <- levels(d$.trt)
    if (! trt_ref %in% lvls_trt)
      abort(sprintf("`trt_ref` does not match a treatment in the data.\nSuitable values are: %s",
                    ifelse(length(lvls_trt) <= 5,
                           paste0(lvls_trt, collapse = ", "),
                           paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))
    d$.trt <- forcats::fct_relevel(d$.trt, trt_ref)
  }

  if (!is.null(.trtclass)) {
    d <- tibble::add_column(d, .trtclass = nfactor(.trtclass))
    class_lookup <- d %>%
      dplyr::distinct(.data$.trt, .data$.trtclass) %>%
      dplyr::arrange(.data$.trt)
    class_ref <- as.character(class_lookup[[1, ".trtclass"]])
    d$.trtclass <- forcats::fct_relevel(d$.trtclass, class_ref)
    classes <- forcats::fct_relevel(nfactor(class_lookup$.trtclass), class_ref)
  } else {
    classes <- NULL
  }

  if (o_type == "continuous") {
    d <- tibble::add_column(d, .y = .y)
  } else if (o_type == "binary") {
    d <- tibble::add_column(d, .r = .r)
  } else if (o_type == "rate") {
    d <- tibble::add_column(d, .r = .r, .E = .E)
  } else if (o_type %in% c("ordered", "competing")) {
    # Store internally as a standard matrix (strip class attribute) as a simple
    # work-around for dplyr/vctrs not using s3 inheritance to find matrix methods
    .r <- unclass(.r)

    d <- tibble::add_column(d, .r = .r)
  } else if (o_type == "survival") {
    d <- tibble::add_column(d, .Surv = .Surv)
  }

  drop_reserved <- setdiff(colnames(data), colnames(d))
  d <- dplyr::bind_cols(d, data[, drop_reserved, drop = FALSE])

  # Drop original study and treatment columns
  d <- drop_original(d, data, enquo(study))
  d <- drop_original(d, data, enquo(trt))
  if (!is.null(.trtclass)) d <- drop_original(d, data, enquo(trt_class))

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = d,
         treatments = forcats::fct_unique(d$.trt),
         classes = classes,
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = NA, agd_contrast = NA, ipd = o_type)),
    class = "nma_data")

  # If trt_ref not specified, mark treatments factor as default, calculate
  # current reference trt
  if (is.null(trt_ref)) {
    trt_ref <- get_default_trt_ref(out)
    trt_sort <- order(forcats::fct_relevel(out$treatments, trt_ref))
    out$treatments <- .default(forcats::fct_relevel(out$treatments, trt_ref)[trt_sort])
    out$ipd$.trt <- forcats::fct_relevel(out$ipd$.trt, trt_ref)
    if (!is.null(.trtclass)) {
      class_ref <- as.character(out$classes[trt_sort[1]])
      if (!is.null(trtclass_original_levels))
        class_ref <- c(class_ref,
                       setdiff(intersect(trtclass_original_levels, levels(out$classes)),
                               class_ref))
      out$ipd$.trtclass <- forcats::fct_relevel(out$ipd$.trtclass, class_ref)
      out$classes <- forcats::fct_relevel(out$classes, class_ref)[trt_sort]
    }
  }

  # Add original_levels attributes (if not NULL)
  attr(out$treatments, "original_levels") <- trt_original_levels
  attr(out$studies, "original_levels") <- study_original_levels
  if (!is.null(.trtclass)) attr(out$classes, "original_levels") <- trtclass_original_levels

  return(out)
}


#' Set up arm-based aggregate data
#'
#' Set up a network containing arm-based aggregate data (AgD), such as event
#' counts or mean outcomes on each arm. Multiple data sources may be combined
#' once created using [combine_network()].
#'
#' @template args-data_common
#' @template args-data_y
#' @template args-data_se
#' @template args-data_rE
#' @param n column of `data` specifying Binomial outcome numerator
#' @param sample_size column of `data` giving the sample size in each arm.
#'   Optional, see details.
#'
#' @return An object of class [nma_data]
#' @export

#' @template args-details_trt_ref
#' @template args-details_sample_size
#' @details
#' If a Binomial outcome is specified and `sample_size` is omitted, `n` will be
#' used as the sample size by default. If a Multinomial outcome is specified and
#' `sample_size` is omitted, the sample size will be determined automatically
#' from the supplied counts by default.
#'
#' @template args-details_mutate
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_contrast()] for
#'   contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @template seealso_nma_data
#' @template ex_smoking_network
#' @examples
#'
#' # Plot network
#' plot(smk_net)
set_agd_arm <- function(data,
                        study,
                        trt,
                        y = NULL, se = NULL,
                        r = NULL, n = NULL, E = NULL,
                        sample_size = NULL,
                        trt_ref = NULL,
                        trt_class = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             classes = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- pull_non_null(data, enquo(study))
  if (is.null(.study)) abort("`study` cannot be NULL")
  check_study(.study)

  if (is.factor(.study)) {
    study_original_levels <- levels(.study)
    .study <- forcats::fct_drop(.study)
  } else {
    study_original_levels <- NULL
  }

  if (missing(trt)) abort("Specify `trt`")
  .trt <- pull_non_null(data, enquo(trt))
  if (is.null(.trt)) abort("`trt` cannot be NULL")
  check_trt(.trt)

  if (is.factor(.trt)) {
    trt_original_levels <- levels(.trt)
    .trt <- forcats::fct_drop(.trt)
  } else {
    trt_original_levels <- NULL
  }

  # Check for single-arm studies
  single_arm_studies <- tibble::tibble(.study, .trt) %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::pull(.data$.study)

  if (length(single_arm_studies)) {
    abort(glue::glue("Single-arm studies are not supported: issue with stud{if (length(single_arm_studies) > 1) 'ies' else 'y'} ",
                     glue::glue_collapse(glue::double_quote(single_arm_studies), sep = ", ", last = " and "), "."))
  }

  # Treatment classes
  .trtclass <- pull_non_null(data, enquo(trt_class))
  if (!is.null(.trtclass)) {
    check_trt_class(.trtclass, .trt)

    if (is.factor(.trtclass)) {
      trtclass_original_levels <- levels(.trtclass)
      .trtclass <- forcats::fct_drop(.trtclass)
    } else {
      trtclass_original_levels <- NULL
    }
  }

  if (!is.null(trt_ref) && length(trt_ref) > 1) abort("`trt_ref` must be length 1.")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .se <- pull_non_null(data, enquo(se))
  .r <- pull_non_null(data, enquo(r))
  .n <- pull_non_null(data, enquo(n))
  .E <- pull_non_null(data, enquo(E))

  check_outcome_continuous(.y, .se, with_se = TRUE)

  if (!is.null(.r) && inherits(.r, c("multi_ordered", "multi_competing"))) {
    # Checks are carried out by multi()
    if (inherits(.r, "multi_competing")) abort("Competing multinomial outcomes are not yet supported.")
  } else {
    check_outcome_count(.r, .n, .E)
  }

  o_type <- get_outcome_type(y = .y, se = .se,
                             r = .r, n = .n, E = .E,
                             Surv = NULL)

  # Pull and check sample size
  .sample_size <- pull_non_null(data, enquo(sample_size))
  if (!is.null(.sample_size)) {
    check_sample_size(.sample_size)
  } else if (o_type == "count") {
    .sample_size <- .n
  } else if (o_type %in% c("ordered", "competing")) {
    # Always stored as exclusive counts
    .sample_size <- rowSums(.r, na.rm = TRUE)
  }
  else inform("Note: Optional argument `sample_size` not provided, some features may not be available (see ?set_agd_arm).")

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt)
  )

  if (!is.null(trt_ref)) {
    trt_ref <- as.character(trt_ref)
    lvls_trt <- levels(d$.trt)
    if (! trt_ref %in% lvls_trt)
      abort(sprintf("`trt_ref` does not match a treatment in the data.\nSuitable values are: %s",
                    ifelse(length(lvls_trt) <= 5,
                           paste0(lvls_trt, collapse = ", "),
                           paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))
    d$.trt <- forcats::fct_relevel(d$.trt, trt_ref)
  }

  if (!is.null(.trtclass)) {
    d <- tibble::add_column(d, .trtclass = nfactor(.trtclass))
    class_lookup <- d %>%
      dplyr::distinct(.data$.trt, .data$.trtclass) %>%
      dplyr::arrange(.data$.trt)
    class_ref <- as.character(class_lookup[[1, ".trtclass"]])
    d$.trtclass <- forcats::fct_relevel(d$.trtclass, class_ref)
    classes <- forcats::fct_relevel(nfactor(class_lookup$.trtclass), class_ref)
  } else {
    classes <- NULL
  }

  if (o_type == "continuous") {
    d <- tibble::add_column(d, .y = .y, .se = .se)
  } else if (o_type == "count") {
    d <- tibble::add_column(d, .r = .r, .n = .n)
  } else if (o_type == "rate") {
    d <- tibble::add_column(d, .r = .r, .E = .E)
  } else if (o_type %in% c("ordered", "competing")) {
    # Store internally as a standard matrix (strip class attribute) as a simple
    # work-around for dplyr/vctrs not using s3 inheritance to find matrix methods
    .r <- unclass(.r)

    d <- tibble::add_column(d, .r = .r)
  }

  if (!is.null(.sample_size)) d <- tibble::add_column(d, .sample_size = .sample_size)

  # Bind in original data
  drop_reserved <- setdiff(colnames(data), colnames(d))
  d <- dplyr::bind_cols(d, data[, drop_reserved, drop = FALSE])

  # Drop original study and treatment columns
  d <- drop_original(d, data, enquo(study))
  d <- drop_original(d, data, enquo(trt))
  if (!is.null(.trtclass)) d <- drop_original(d, data, enquo(trt_class))

  # Produce nma_data object
  out <- structure(
    list(agd_arm = d,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         classes = classes,
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = o_type, agd_contrast = NA, ipd = NA)),
    class = "nma_data")

  # If trt_ref not specified, mark treatments factor as default, calculate
  # current reference trt
  if (is.null(trt_ref)) {
    trt_ref <- get_default_trt_ref(out)
    trt_sort <- order(forcats::fct_relevel(out$treatments, trt_ref))
    out$treatments <- .default(forcats::fct_relevel(out$treatments, trt_ref)[trt_sort])
    out$agd_arm$.trt <- forcats::fct_relevel(out$agd_arm$.trt, trt_ref)
    if (!is.null(.trtclass)) {
      class_ref <- as.character(out$classes[trt_sort[1]])
      if (!is.null(trtclass_original_levels))
        class_ref <- c(class_ref,
                       setdiff(intersect(trtclass_original_levels, levels(out$classes)),
                               class_ref))
      out$agd_arm$.trtclass <- forcats::fct_relevel(out$agd_arm$.trtclass, class_ref)
      out$classes <- forcats::fct_relevel(out$classes, class_ref)[trt_sort]
    }
  }

  # Add original_levels attributes (if not NULL)
  attr(out$treatments, "original_levels") <- trt_original_levels
  attr(out$studies, "original_levels") <- study_original_levels
  if (!is.null(.trtclass)) attr(out$classes, "original_levels") <- trtclass_original_levels

  return(out)
}


#' Set up contrast-based aggregate data
#'
#' Set up a network containing contrast-based aggregate data (AgD), i.e.
#' summaries of relative effects between treatments such as log Odds Ratios.
#' Multiple data sources may be combined once created using [combine_network()].
#'
#' @template args-data_common
#' @template args-data_y
#' @template args-data_se
#' @param sample_size column of `data` giving the sample size in each arm.
#'   Optional, see details.
#'
#' @details Each study should have a single reference/baseline treatment,
#'   against which relative effects in the other arm(s) are given. For the
#'   reference arm, include a data row with continuous outcome `y` equal to
#'   `NA`. If a study has three or more arms (so two or more relative effects),
#'   set the standard error `se` for the reference arm data row equal to the
#'   standard error of the mean outcome on the reference arm (this determines
#'   the covariance of the relative effects, when expressed as differences in
#'   mean outcomes between arms).
#'
#' @template args-details_mutate
#'
#' @template args-details_trt_ref
#' @template args-details_sample_size
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_arm()] for
#'   arm-based aggregate data, and [combine_network()] for combining several
#'   data sources in one network.
#' @template seealso_nma_data
#' @examples
#' # Set up network of Parkinson's contrast data
#' head(parkinsons)
#'
#' park_net <- set_agd_contrast(parkinsons,
#'                              study = studyn,
#'                              trt = trtn,
#'                              y = diff,
#'                              se = se_diff,
#'                              sample_size = n)
#'
#' # Print details
#' park_net
#'
#' # Plot network
#' plot(park_net)

set_agd_contrast <- function(data,
                             study,
                             trt,
                             y = NULL, se = NULL,
                             sample_size = NULL,
                             trt_ref = NULL,
                             trt_class = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             classes = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }


  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- pull_non_null(data, enquo(study))
  if (is.null(.study)) abort("`study` cannot be NULL")
  check_study(.study)

  if (is.factor(.study)) {
    study_original_levels <- levels(.study)
    .study <- forcats::fct_drop(.study)
  } else {
    study_original_levels <- NULL
  }

  if (missing(trt)) abort("Specify `trt`")
  .trt <- pull_non_null(data, enquo(trt))
  if (is.null(.trt)) abort("`trt` cannot be NULL")
  check_trt(.trt)

  if (is.factor(.trt)) {
    trt_original_levels <- levels(.trt)
    .trt <- forcats::fct_drop(.trt)
  } else {
    trt_original_levels <- NULL
  }

  # Check for single-arm studies
  single_arm_studies <- tibble::tibble(.study, .trt) %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::pull(.data$.study)

  if (length(single_arm_studies)) {
    abort(glue::glue("Single-arm studies are not supported: issue with stud{if (length(single_arm_studies) > 1) 'ies' else 'y'} ",
                     glue::glue_collapse(glue::double_quote(single_arm_studies), sep = ", ", last = " and "), "."))
  }

  # Treatment classes
  .trtclass <- pull_non_null(data, enquo(trt_class))
  if (!is.null(.trtclass)) {
    check_trt_class(.trtclass, .trt)

    if (is.factor(.trtclass)) {
      trtclass_original_levels <- levels(.trtclass)
      .trtclass <- forcats::fct_drop(.trtclass)
    } else {
      trtclass_original_levels <- NULL
    }
  }

  if (!is.null(trt_ref) && length(trt_ref) > 1) abort("`trt_ref` must be length 1.")

  # Pull and check outcomes
  .y <- pull_non_null(data, enquo(y))
  .se <- pull_non_null(data, enquo(se))

  if (is.null(.y)) abort("Specify continuous outcome `y`")
  if (rlang::is_list(.y) || !is.null(dim(.y)))
    abort("Continuous outcome `y` must be a regular column (not a list or matrix column)")

  if (is.null(.se)) abort("Specify standard error `se`")
  if (rlang::is_list(.se) || !is.null(dim(.se)))
    abort("Standard error `se` must be a regular column (not a list or matrix column)")

  # Pull and check sample size
  .sample_size <- pull_non_null(data, enquo(sample_size))
  if (!is.null(.sample_size)) {
    check_sample_size(.sample_size)
  } else {
    inform("Note: Optional argument `sample_size` not provided, some features may not be available (see ?set_agd_contrast).")
  }

  # Determine baseline arms by .y = NA
  bl <- is.na(.y)

  # if (anyDuplicated(.study[bl])) abort("Multiple baseline arms (where y = NA) for a study.")

  tibble::tibble(.study, .trt, bl, .se) %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::mutate(n_arms = dplyr::n(),
                  n_bl = sum(.data$bl)) %>%
    {
      if (any(.$n_bl > 1))
        abort("Multiple baseline arms (where y = NA) in a study or studies.")
      else if (any(.$n_bl == 0))
        abort("Study or studies without a specified baseline arm (where y = NA).")
      else .
    } %>%
    dplyr::filter(.data$bl, .data$n_arms > 2) %>%
    {
      check_outcome_continuous(1, .$.se, with_se = TRUE,
                               append = " on baseline arms in studies with >2 arms.")
    }

  check_outcome_continuous(.y[!bl], .se[!bl], with_se = TRUE,
    append = " for non-baseline rows (i.e. those specifying contrasts against baseline).")

  o_type <- get_outcome_type(y = .y[!bl], se = .se[!bl],
                             r = NULL, n = NULL, E = NULL, Surv = NULL)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt),
    .y = .y,
    .se = .se)

  if (!is.null(trt_ref)) {
    trt_ref <- as.character(trt_ref)
    lvls_trt <- levels(d$.trt)
    if (! trt_ref %in% lvls_trt)
      abort(sprintf("`trt_ref` does not match a treatment in the data.\nSuitable values are: %s",
                    ifelse(length(lvls_trt) <= 5,
                           paste0(lvls_trt, collapse = ", "),
                           paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))
    d$.trt <- forcats::fct_relevel(d$.trt, trt_ref)
  }

  if (!is.null(.trtclass)) {
    d <- tibble::add_column(d, .trtclass = nfactor(.trtclass))
    class_lookup <- d %>%
      dplyr::distinct(.data$.trt, .data$.trtclass) %>%
      dplyr::arrange(.data$.trt)
    class_ref <- as.character(class_lookup[[1, ".trtclass"]])
    d$.trtclass <- forcats::fct_relevel(d$.trtclass, class_ref)
    classes <- forcats::fct_relevel(nfactor(class_lookup$.trtclass), class_ref)
  } else {
    classes <- NULL
  }

  if (!is.null(.sample_size)) {
    d <- tibble::add_column(d, .sample_size = .sample_size)
  }

  # Bind in original data
  drop_reserved <- setdiff(colnames(data), colnames(d))
  d <- dplyr::bind_cols(d, data[, drop_reserved, drop = FALSE])

  # Drop original study and treatment columns
  d <- drop_original(d, data, enquo(study))
  d <- drop_original(d, data, enquo(trt))
  if (!is.null(.trtclass)) d <- drop_original(d, data, enquo(trt_class))

  # Make sure rows from each study are next to each other (required for Stan resdev/log_lik code)
  d <- dplyr::mutate(d, .study_inorder = forcats::fct_inorder(.data$.study)) %>%
    dplyr::arrange(.data$.study_inorder) %>%
    dplyr::select(-".study_inorder")

  # Check covariance matrices are positive definite
  agd_Sigma <- make_Sigma(d)
  posdef <- purrr::map_lgl(agd_Sigma, ~all(eigen(., only.values = TRUE)$values > sqrt(.Machine$double.eps)))

  if (any(!posdef)) {
    abort(glue::glue("Contrast covariance matrix not positive definite for stud{if (sum(!posdef) > 1) 'ies' else 'y'} ",
                     glue::glue_collapse(glue::double_quote(names(posdef[!posdef])), sep = ", ", last = " and "), ".\n",
                     "Check the `se` column."))
  }

  # Produce nma_data object
  out <- structure(
    list(agd_arm = NULL,
         agd_contrast = d,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         classes = classes,
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = NA, agd_contrast = o_type, ipd = NA)),
    class = "nma_data")

  # If trt_ref not specified, mark treatments factor as default, calculate
  # current reference trt
  if (is.null(trt_ref)) {
    trt_ref <- get_default_trt_ref(out)
    trt_sort <- order(forcats::fct_relevel(out$treatments, trt_ref))
    out$treatments <- .default(forcats::fct_relevel(out$treatments, trt_ref)[trt_sort])
    out$agd_contrast$.trt <- forcats::fct_relevel(out$agd_contrast$.trt, trt_ref)
    if (!is.null(.trtclass)) {
      class_ref <- as.character(out$classes[trt_sort[1]])
      if (!is.null(trtclass_original_levels))
        class_ref <- c(class_ref,
                       setdiff(intersect(trtclass_original_levels, levels(out$classes)),
                               class_ref))
      out$agd_contrast$.trtclass <- forcats::fct_relevel(out$agd_contrast$.trtclass, class_ref)
      out$classes <- forcats::fct_relevel(out$classes, class_ref)[trt_sort]
    }
  }

  # Add original_levels attributes (if not NULL)
  attr(out$treatments, "original_levels") <- trt_original_levels
  attr(out$studies, "original_levels") <- study_original_levels
  if (!is.null(.trtclass)) attr(out$classes, "original_levels") <- trtclass_original_levels

  return(out)
}


#' Set up aggregate survival data
#'
#' Set up a network containing aggregate survival data (AgD) in the form of
#' event/censoring times (e.g. reconstructed from digitized Kaplan-Meier curves)
#' and covariate summary statistics from each study. Multiple data sources may be
#' combined once created using [combine_network()].
#'
#' @template args-data_common
#' @template args-data_Surv
#' @param covariates data frame of covariate summary statistics for each study
#'   or study arm, with corresponding `study` and `trt` columns to match to
#'   those in `data`
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @template args-details_trt_ref
#' @template args-details_mutate
#'
#' @seealso [set_ipd()] for individual patient data, [set_agd_contrast()] for
#'   contrast-based aggregate data, and [combine_network()] for combining
#'   several data sources in one network.
#' @template seealso_nma_data
#' @examples
#' ## Newly diagnosed multiple myeloma
#'
#' head(ndmm_agd)  # Reconstructed Kaplan-Meier data
#' ndmm_agd_covs   # Summary covariate information on each arm
#'
#' set_agd_surv(ndmm_agd,
#'              study = studyf,
#'              trt = trtf,
#'              Surv = Surv(eventtime, status),
#'              covariates = ndmm_agd_covs)
#'
set_agd_surv <- function(data,
                         study,
                         trt,
                         Surv,
                         covariates = NULL,
                         trt_ref = NULL,
                         trt_class = NULL) {

  # Check data is data frame
  if (!inherits(data, "data.frame")) abort("Argument `data` should be a data frame")
  if (nrow(data) == 0) {
    return(
      structure(
        list(agd_arm = NULL,
             agd_contrast = NULL,
             ipd = NULL,
             treatments = NULL,
             classes = NULL,
             studies = NULL),
        class = "nma_data")
    )
  }

  if (!is.null(covariates)) {
    if (!inherits(covariates, "data.frame")) abort("Argument `covariates` should be a data frame")
  }

  # Pull study and treatment columns
  if (missing(study)) abort("Specify `study`")
  .study <- pull_non_null(data, enquo(study))
  if (is.null(.study)) abort("`study` cannot be NULL")
  check_study(.study)

  if (is.factor(.study)) {
    study_original_levels <- levels(.study)
    .study <- forcats::fct_drop(.study)
  } else {
    study_original_levels <- NULL
  }

  if (missing(trt)) abort("Specify `trt`")
  .trt <- pull_non_null(data, enquo(trt))
  if (is.null(.trt)) abort("`trt` cannot be NULL")
  check_trt(.trt)

  if (is.factor(.trt)) {
    trt_original_levels <- levels(.trt)
    .trt <- forcats::fct_drop(.trt)
  } else {
    trt_original_levels <- NULL
  }

  # Check for single-arm studies
  single_arm_studies <- tibble::tibble(.study, .trt) %>%
    dplyr::distinct(.data$.study, .data$.trt) %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::pull(.data$.study)

  if (length(single_arm_studies)) {
    inform(glue::glue("Single-arm stud{if (length(single_arm_studies) > 1) 'ies' else 'y'} present in the network: ",
                      glue::glue_collapse(glue::double_quote(as.character(single_arm_studies)), sep = ", ", last = " and "), "."))
  }

  # Treatment classes
  .trtclass <- pull_non_null(data, enquo(trt_class))
  if (!is.null(.trtclass)) {
    check_trt_class(.trtclass, .trt)

    if (is.factor(.trtclass)) {
      trtclass_original_levels <- levels(.trtclass)
      .trtclass <- forcats::fct_drop(.trtclass)
    } else {
      trtclass_original_levels <- NULL
    }
  }

  if (!is.null(trt_ref) && length(trt_ref) > 1) abort("`trt_ref` must be length 1.")

  # Pull and check outcomes
  .Surv <- pull_non_null(data, enquo(Surv))

  check_outcome_survival(.Surv)

  o_type <- get_outcome_type(y = NULL, se = NULL,
                             r = NULL, n = NULL, E = NULL,
                             Surv = .Surv)

  # Create tibble in standard format
  d <- tibble::tibble(
    .study = nfactor(.study),
    .trt = nfactor(.trt)
  )

  if (!is.null(trt_ref)) {
    trt_ref <- as.character(trt_ref)
    lvls_trt <- levels(d$.trt)
    if (! trt_ref %in% lvls_trt)
      abort(sprintf("`trt_ref` does not match a treatment in the data.\nSuitable values are: %s",
                    ifelse(length(lvls_trt) <= 5,
                           paste0(lvls_trt, collapse = ", "),
                           paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))
    d$.trt <- forcats::fct_relevel(d$.trt, trt_ref)
  }

  if (!is.null(.trtclass)) {
    d <- tibble::add_column(d, .trtclass = nfactor(.trtclass))
    class_lookup <- d %>%
      dplyr::distinct(.data$.trt, .data$.trtclass) %>%
      dplyr::arrange(.data$.trt)
    class_ref <- as.character(class_lookup[[1, ".trtclass"]])
    d$.trtclass <- forcats::fct_relevel(d$.trtclass, class_ref)
    classes <- forcats::fct_relevel(nfactor(class_lookup$.trtclass), class_ref)
  } else {
    classes <- NULL
  }

  d <- tibble::add_column(d, .Surv = .Surv)

  # Add in sample size
  d <- dplyr::group_by(d, .data$.study, .data$.trt) %>%
    dplyr::mutate(.sample_size = dplyr::n()) %>%
    dplyr::ungroup()

  # Bind in original data
  drop_reserved <- setdiff(colnames(data), colnames(d))
  d <- dplyr::bind_cols(d, data[, drop_reserved, drop = FALSE])

  # Drop original study and treatment columns
  d <- drop_original(d, data, enquo(study))
  d <- drop_original(d, data, enquo(trt))
  if (!is.null(.trtclass)) d <- drop_original(d, data, enquo(trt_class))

  # Nest survival data, keeping original data nested for later too
  d <- dplyr::group_by(d, dplyr::across(dplyr::any_of(c(".study", ".trt", ".trtclass", ".sample_size")))) %>%
    tidyr::nest(.Surv = ".Surv",
                .data_orig = !dplyr::any_of(c(".study", ".trt", ".trtclass", ".sample_size", ".Surv"))) %>%
    dplyr::ungroup()

  # Join covariate details
  if (!is.null(covariates)) {
    covariates <- dplyr::ungroup(covariates)

    .cov_study <- pull_non_null(covariates, enquo(study))
    if (is.null(.cov_study)) abort("`study` cannot be NULL")
    check_study(.cov_study)
    .cov_study <- forcats::fct_drop(nfactor(.cov_study))

    .cov_trt <- pull_non_null(covariates, enquo(trt))
    if (is.null(.cov_trt)) abort("`trt` cannot be NULL")
    check_trt(.cov_trt)
    .cov_trt <- forcats::fct_drop(nfactor(.cov_trt))

    if (!is.null(trt_ref)) {
      .cov_trt <- forcats::fct_relevel(.cov_trt, trt_ref)
    }

    covs <- dplyr::mutate(covariates, .study = .cov_study, .trt = .cov_trt)
    covs <- drop_original(covs, data, enquo(study))
    covs <- drop_original(covs, data, enquo(trt))

    d <- tryCatch(
          dplyr::inner_join(d, covs, by = c(".study", ".trt"),
                            multiple = "first",
                            unmatched = c("error", "drop")),
          error = function(e) abort("Not all study arms in `data` have matching rows in `covariates`", parent = e))

    # Re-drop factors, in case extra unneeded trt/study rows included in covariate data
    d$.trt <- forcats::fct_drop(d$.trt)
    d$.study <- forcats::fct_drop(d$.study)
  }

  # Produce nma_data object
  out <- structure(
    list(agd_arm = d,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = forcats::fct_unique(d$.trt),
         classes = classes,
         studies = forcats::fct_unique(d$.study),
         outcome = list(agd_arm = o_type, agd_contrast = NA, ipd = NA)),
    class = "nma_data")

  # If trt_ref not specified, mark treatments factor as default, calculate
  # current reference trt
  if (is.null(trt_ref)) {
    trt_ref <- get_default_trt_ref(out)
    trt_sort <- order(forcats::fct_relevel(out$treatments, trt_ref))
    out$treatments <- .default(forcats::fct_relevel(out$treatments, trt_ref)[trt_sort])
    out$agd_arm$.trt <- forcats::fct_relevel(out$agd_arm$.trt, trt_ref)
    if (!is.null(.trtclass)) {
      class_ref <- as.character(out$classes[trt_sort[1]])
      if (!is.null(trtclass_original_levels))
        class_ref <- c(class_ref,
                       setdiff(intersect(trtclass_original_levels, levels(out$classes)),
                               class_ref))
      out$agd_arm$.trtclass <- forcats::fct_relevel(out$agd_arm$.trtclass, class_ref)
      out$classes <- forcats::fct_relevel(out$classes, class_ref)[trt_sort]
    }
  }

  # Add original_levels attributes (if not NULL)
  attr(out$treatments, "original_levels") <- trt_original_levels
  attr(out$studies, "original_levels") <- study_original_levels
  if (!is.null(.trtclass)) attr(out$classes, "original_levels") <- trtclass_original_levels

  return(out)

}


#' Combine multiple data sources into one network
#'
#' Multiple data sources created using [set_ipd()], [set_agd_arm()], or
#' [set_agd_contrast()] can be combined into a single network for analysis.
#'
#' @param ... multiple data sources, as defined using the `set_*` functions
#' @param trt_ref reference treatment for the entire network, as a string (or
#'   coerced as such) referring to the levels of the treatment factor variable
#'
#' @return An object of class [nma_data]
#' @export
#'
#' @seealso [set_ipd()], [set_agd_arm()], and [set_agd_contrast()] for defining
#'   different data sources.
#' @template seealso_nma_data
#'
#' @examples ## Parkinson's - combining contrast- and arm-based data
#' studies <- parkinsons$studyn
#' (parkinsons_arm <- parkinsons[studies %in% 1:3, ])
#' (parkinsons_contr <- parkinsons[studies %in% 4:7, ])
#'
#' park_arm_net <- set_agd_arm(parkinsons_arm,
#'                             study = studyn,
#'                             trt = trtn,
#'                             y = y,
#'                             se = se,
#'                             sample_size = n)
#'
#' park_contr_net <- set_agd_contrast(parkinsons_contr,
#'                                    study = studyn,
#'                                    trt = trtn,
#'                                    y = diff,
#'                                    se = se_diff,
#'                                    sample_size = n)
#'
#' park_net <- combine_network(park_arm_net, park_contr_net)
#'
#' # Print network details
#' park_net
#'
#' # Plot network
#' plot(park_net, weight_edges = TRUE, weight_nodes = TRUE)
#'
#' @examples ## Plaque Psoriasis - combining IPD and AgD in a network
#' @template ex_plaque_psoriasis_network
#' @examples
#'
#' # Plot network
#' plot(pso_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE)
combine_network <- function(..., trt_ref) {
  s <- list(...)

  # Check that arguments all inherit from nma_data class
  if (!purrr::every(s, inherits, what = "nma_data")) {
    abort("Expecting to combine objects of class `nma_data`, created using set_* functions")
  }

  # Combine treatment code factor
  trts <- stringr::str_sort(forcats::lvls_union(purrr::map(s, "treatments")), numeric = TRUE)

  # Use original factor order (possibly with more levels) if available
  trt_original_levels <- purrr::map(purrr::map(s, "treatments"), attr, which = "original_levels")
  if (!any(purrr::map_lgl(trt_original_levels, is.null)) &&
      all(purrr::map_lgl(trt_original_levels, ~identical(., trt_original_levels[[1]])))) {

    trt_original_levels <- trt_original_levels[[1]]
    trts <- intersect(trt_original_levels, trts)
  } else {
    trt_original_levels <- NULL
  }

  if (!missing(trt_ref)) {
    if (! trt_ref %in% trts) {
      abort(sprintf("`trt_ref` does not match a treatment in the network.\nSuitable values are: %s",
                      ifelse(length(trts) <= 5,
                             paste0(trts, collapse = ", "),
                             paste0(paste0(trts[1:5], collapse = ", "), ", ..."))))
    }
    trts <- c(trt_ref, setdiff(trts, trt_ref))
  }

  # Combine classes factor
  has_classes <- purrr::map_lgl(purrr::map(s, "classes"), ~!is.null(.))

  if (all(has_classes)) {
    class_lookup <- tibble::tibble(.trt = forcats::fct_c(!!! purrr::map(s, "treatments")),
                                   .trtclass = forcats::fct_c(!!! purrr::map(s, "classes"))) %>%
      dplyr::mutate(.trt = forcats::fct_relevel(.data$.trt, trts)) %>%
      dplyr::distinct(.data$.trt, .data$.trtclass) %>%
      dplyr::arrange(.data$.trt)

    check_trt_class(class_lookup$.trtclass, class_lookup$.trt)

    class_lvls <- stringr::str_sort(levels(class_lookup$.trtclass), numeric = TRUE)

    # Use original factor order (possibly with more levels) if available
    class_original_levels <- purrr::map(purrr::map(s, "classes"), attr, which = "original_levels")
    if (!any(purrr::map_lgl(class_original_levels, is.null)) &&
        all(purrr::map_lgl(class_original_levels, ~identical(., class_original_levels[[1]])))) {

      class_original_levels <- class_original_levels[[1]]
      class_lvls <- intersect(class_original_levels, class_lvls)
    } else {
      class_original_levels <- NULL
    }

    class_ref <- as.character(class_lookup[[1, ".trtclass"]])
    class_lvls <- c(class_ref, setdiff(class_lvls, class_ref))

    class_lookup$.trtclass <- forcats::fct_relevel(class_lookup$.trtclass, class_lvls)

    classes <- class_lookup$.trtclass
  } else if (any(has_classes)) {
    warn("Not all data sources have defined treatment classes. Removing treatment class information.")
    classes <- NULL
  } else {
    classes <- NULL
  }

  # Check that no studies are duplicated between data sources
  all_studs <- purrr::flatten_chr(purrr::map(s, ~levels(.$studies)))
  if (anyDuplicated(all_studs)) {
    abort(sprintf("Studies with same label found in multiple data sources: %s",
                  paste0(unique(all_studs[duplicated(all_studs)]), collapse = ", ")))
  }

  # Combine study code factor
  studs <- stringr::str_sort(forcats::lvls_union(purrr::map(s, "studies")), numeric = TRUE)

  # Use original factor order (possibly with more levels) if available
  study_original_levels <- purrr::map(purrr::map(s, "studies"), attr, which = "original_levels")
  if (!any(purrr::map_lgl(study_original_levels, is.null)) &&
      all(purrr::map_lgl(study_original_levels, ~identical(., study_original_levels[[1]])))) {

    study_original_levels <- study_original_levels[[1]]
    studs <- intersect(study_original_levels, studs)
  } else {
    study_original_levels <- NULL
  }

  # Get ipd
  ipd <- purrr::map(s, "ipd")
  if (!rlang::is_empty(ipd)) {
    for (j in 1:length(ipd)) {
      if (rlang::is_empty(ipd[[j]])) next
      ipd[[j]]$.trt <- forcats::lvls_expand(ipd[[j]]$.trt, trts)
      ipd[[j]]$.study <- forcats::lvls_expand(ipd[[j]]$.study, studs)
      if (!is.null(classes)) ipd[[j]]$.trtclass <- forcats::lvls_expand(ipd[[j]]$.trtclass, class_lvls)
    }
  }

  # Get agd_arm
  agd_arm <- purrr::map(s, "agd_arm")
  if (!rlang::is_empty(agd_arm)) {
    for (j in 1:length(agd_arm)) {
      if (rlang::is_empty(agd_arm[[j]])) next
      agd_arm[[j]]$.trt <- forcats::lvls_expand(agd_arm[[j]]$.trt, trts)
      agd_arm[[j]]$.study <- forcats::lvls_expand(agd_arm[[j]]$.study, studs)
      if (!is.null(classes))
        agd_arm[[j]]$.trtclass <- forcats::lvls_expand(agd_arm[[j]]$.trtclass, class_lvls)
    }
  }

  # Get agd_contrast
  agd_contrast <- purrr::map(s, "agd_contrast")
  if (!rlang::is_empty(agd_contrast)) {
    for (j in 1:length(agd_contrast)) {
      if (rlang::is_empty(agd_contrast[[j]])) next
      agd_contrast[[j]]$.trt <- forcats::lvls_expand(agd_contrast[[j]]$.trt, trts)
      agd_contrast[[j]]$.study <- forcats::lvls_expand(agd_contrast[[j]]$.study, studs)
      if (!is.null(classes))
        agd_contrast[[j]]$.trtclass <- forcats::lvls_expand(agd_contrast[[j]]$.trtclass, class_lvls)
    }
  }

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

  # Additional checks when combining multinomial outcomes
  if (o_ipd %in% c("ordered", "competing"))
    check_multi_combine(purrr::map(ipd, ".r"))
  if (o_agd_arm %in% c("ordered", "competing"))
    check_multi_combine(purrr::map(agd_arm, ".r"))
  if (o_ipd %in% c("ordered", "competing") && o_agd_arm %in% c("ordered", "competing"))
    check_multi_combine(purrr::map(c(ipd, agd_arm), ".r"))

  # Check integration points
  if (any(purrr::map_lgl(s, inherits, what = "mlnmr_data")) &&
      all(purrr::map_lgl(s, inherits, what = "mlnmr_data") |
        (purrr::map_lgl(s, has_ipd) & !purrr::map_lgl(s, has_agd_arm) & !purrr::map_lgl(s, has_agd_contrast)))) {

    has_int <- TRUE

    s_int <- s[purrr::map_lgl(s, inherits, what = "mlnmr_data")]

    # Check int_names
    l_int_names <- purrr::map(s_int, "int_names")
    if (!all(purrr::map_lgl(l_int_names, ~identical(., l_int_names[[1]])))) {
      abort("Cannot combine AgD sources with different sets of covariates with integration points.")
    } else {
      int_names <- l_int_names[[1]]
    }

    # Check n_int
    l_n_int <- purrr::map_dbl(s_int, "n_int")
    if (!all(l_n_int == l_n_int[1])) {
      abort("Cannot combine AgD sources with different numbers of integration points.")
    } else {
      n_int <- l_n_int[[1]]
    }

  } else if (any(purrr::map_lgl(s, inherits, what = "mlnmr_data"))) {
    message("Dropping integration points from combined network, some AgD sources do not have integration points.")
    has_int <- FALSE
  } else {
    has_int <- FALSE
  }

  # Produce nma_data object
  ipd <- dplyr::bind_rows(ipd)
  agd_arm <- dplyr::bind_rows(agd_arm)
  agd_contrast <- dplyr::bind_rows(agd_contrast)

  out <- structure(
    list(agd_arm = agd_arm,
         agd_contrast = agd_contrast,
         ipd = ipd,
         treatments = factor(trts, levels = trts),
         classes = classes,
         studies = factor(studs, levels = studs),
         outcome = outcome),
    class = "nma_data")

  # Integration data setup
  if (has_int) {
    out$n_int <- n_int
    out$int_names <- int_names

    class(out) <- c("mlnmr_data", class(out))
  }

  # If trt_ref not specified, mark treatments factor as default, calculate
  # current reference trt
  if (missing(trt_ref)) {
    trt_ref <- get_default_trt_ref(out)
    trt_sort <- order(forcats::fct_relevel(out$treatments, trt_ref))
    out$treatments <- .default(forcats::fct_relevel(out$treatments, trt_ref)[trt_sort])

    if (has_ipd(out))
      out$ipd$.trt <- forcats::fct_relevel(out$ipd$.trt, trt_ref)
    if (has_agd_arm(out))
      out$agd_arm$.trt <- forcats::fct_relevel(out$agd_arm$.trt, trt_ref)
    if (has_agd_contrast(out))
      out$agd_contrast$.trt <- forcats::fct_relevel(out$agd_contrast$.trt, trt_ref)

    if (!is.null(classes)) {
      class_ref <- as.character(out$classes[trt_sort[1]])
      out$classes <- forcats::fct_relevel(out$classes, class_ref)[trt_sort]

      if (has_ipd(out))
        out$ipd$.trtclass <- forcats::fct_relevel(out$ipd$.trtclass, class_ref)
      if (has_agd_arm(out))
        out$agd_arm$.trtclass <- forcats::fct_relevel(out$agd_arm$.trtclass, class_ref)
      if (has_agd_contrast(out))
        out$agd_contrast$.trtclass <- forcats::fct_relevel(out$agd_contrast$.trtclass, class_ref)
    }
  }

  # Add original_levels attributes (if not NULL)
  attr(out$treatments, "original_levels") <- trt_original_levels
  attr(out$studies, "original_levels") <- study_original_levels
  if (!is.null(classes)) attr(out$classes, "original_levels") <- class_original_levels

  return(out)
}

#' Multinomial outcome data
#'
#' This function aids the specification of multinomial outcome data when setting
#' up a network with [set_agd_arm()] or [set_ipd()]. It takes a set of columns
#' (or, more generally, numeric vectors of the same length) of outcome counts in
#' each category, and binds these together to produce a matrix.
#'
#' @param ... Two or more numeric columns (or vectors) of category counts.
#'   Argument names (optional) will be used to label the categories.
#' @param inclusive Logical, are ordered category counts inclusive (`TRUE`) or
#'   exclusive (`FALSE`)? Default `FALSE`. Only used when `type = "ordered"`.
#'   See details.
#' @param type String, indicating whether categories are `"ordered"` or
#'   `"competing"`. Currently only ordered categorical outcomes are supported by
#'   the modelling functions in this package.
#'
#' @details When specifying ordered categorical counts, these can either be
#'   given as *exclusive* counts (`inclusive = FALSE`, the default) where
#'   individuals are only counted in the highest category they achieve, or
#'   *inclusive* counts (`inclusive = TRUE`) where individuals are counted in
#'   every category up to and including the highest category achieved.
#'   (Competing outcomes, by nature, are always specified as exclusive counts.)
#'
#'   `NA` values can be used to indicate categories/cutpoints that were not
#'   measured.
#'
#' @return A matrix of (exclusive) category counts
#' @export
#'
#' @examples
#' # These two data sets specify the same ordered categorical data for outcomes
#' # r0 < r1 < r2, but the first uses the "inclusive" format and the second the
#' # "exclusive" format.
#' df_inclusive <- tibble::tribble(~r0, ~r1, ~r2,
#'                                 1, 1, 1,
#'                                 5, 4, 1,
#'                                 5, 2, 2,
#'                                 10, 5, 0,
#'                                 5, 5, 0,
#'                                 7, NA, 6,   # Achieved r2 or not (no r1)
#'                                 10, 4, NA)  # Achieved r1 or not (no r2)
#'
#' df_exclusive <- tibble::tribble(~r0, ~r1, ~r2,
#'                                 0, 0, 1,
#'                                 1, 3, 1,
#'                                 3, 0, 2,
#'                                 5, 5, 0,
#'                                 0, 5, 0,
#'                                 1, NA, 6,   # Achieved r2 or not (no r1)
#'                                 6, 4, NA)   # Achieved r1 or not (no r2)
#'
#' (r_inclusive <- with(df_inclusive, multi(r0, r1, r2, inclusive = TRUE)))
#' (r_exclusive <- with(df_exclusive, multi(r0, r1, r2, inclusive = FALSE)))
#'
#' # Counts are always stored in exclusive format
#' stopifnot(isTRUE(all.equal(r_inclusive, r_exclusive)))
#'
#'
#' ## HTA Plaque Psoriasis
#' library(dplyr)
#'
#' # Ordered outcomes here are given as "exclusive" counts
#' head(hta_psoriasis)
#'
#' # Calculate lowest category count (failure to achieve PASI 50)
#' pso_dat <- hta_psoriasis %>%
#'   mutate(`PASI<50` = sample_size - rowSums(cbind(PASI50, PASI75, PASI90), na.rm = TRUE))
#'
#' # Set up network
#' pso_net <- set_agd_arm(pso_dat,
#'                        study = paste(studyc, year),
#'                        trt = trtc,
#'                        r = multi(`PASI<50`, PASI50, PASI75, PASI90,
#'                                  inclusive = FALSE,
#'                                  type = "ordered"))
#'
#' pso_net
#'
multi <- function(..., inclusive = FALSE, type = c("ordered", "competing")) {
  if (packageVersion("dplyr") < "1.0.0")
    abort("Multinomial outcomes require `dplyr` package version 1.0.0 or later.")

  # Argument checks
  if (!rlang::is_bool(inclusive)) abort("`inclusive` must be a logical value TRUE/FALSE")
  type <- rlang::arg_match(type)

  if (type == "competing" && inclusive) {
    warn("Ignoring inclusive = TRUE, competing outcomes are always given by exclusive counts.")
    inclusive <- FALSE
  }

  # Collect dots
  # We take this route via quosures to simplify automatic naming
  q_dots <- rlang::enquos(..., .named = TRUE)

  if (length(q_dots) < 2) abort("At least 2 outcomes must be specified in `...`")

  if (anyDuplicated(names(q_dots))) {
    dups <- unique(names(q_dots)[duplicated(names(q_dots))])
    abort(glue::glue("Duplicate outcome category labels ",
                     glue::glue_collapse(glue::double_quote(dups), sep = ", ", last = " and "),
                     "."))
  }

  # Construct outcome matrix
  dots <- purrr::map(q_dots, rlang::eval_tidy)

  # Input vectors must be same length (or length 1 for recycling)
  dots_lengths <- lengths(dots)
  if (length(unique(dots_lengths[dots_lengths > 1])) > 1) abort("Input vectors in `...` must be the same length (or length 1).")

  out <- do.call(cbind, dots)

  # Check counts
  if (!is.numeric(out)) abort("Categorical outcome count must be numeric")
  if (any(is.nan(out))) abort("Categorical outcome count cannot be NaN")
  if (any(is.infinite(out))) abort("Categorical outcome count cannot be Inf")
  if (!rlang::is_integerish(out)) abort("Categorical outcome count must be integer-valued")
  if (any(out < 0, na.rm = TRUE)) abort("Categorical outcome count must be non-negative")

  if (type == "ordered") {
    if (any(c1_na <- is.na(out[, 1]))) {
      abort(glue::glue("Ordered outcome counts cannot be missing in the lowest category.\n",
                       "NAs found in row{if (sum(c1_na) > 1) 's' else ''} ",
                       glue::glue_collapse(which(c1_na), sep = ", ", width = 30, last = " and "), "."))
    }
  }


  if (any(only1 <- apply(out, 1, function(x) sum(!is.na(x))) < 2)) {
    abort(glue::glue("Outcome counts must be present for at least 2 categories.\n",
                     "Issues in row{if (sum(only1) > 1) 's' else ''} ",
                     glue::glue_collapse(which(only1), sep = ", ", width = 30, last = " and "), "."))
  }

  if (inclusive) {
    if (any(non_decreasing <- apply(out, 1, function(x) max(diff(x[!is.na(x)]))) > 0)) {
      abort(glue::glue("Inclusive ordered outcome counts must be decreasing or constant across increasing categories.\n",
                       "Increasing counts found in row{if (sum(non_decreasing) > 1) 's' else ''} ",
                       glue::glue_collapse(which(non_decreasing), sep = ", ", width = 30, last = " and "), "."))
    }

    # Store internally as exclusive counts
    ncat <- ncol(out)
    for (i in 1:nrow(out)) {
      j <- 1L
      k <- 2L
      while (k <= ncat) {
        if (is.na(out[i, k])) {
          k <- k + 1L
        } else {
          out[i, j] <- out[i, j] - out[i, k]
          j <- k
          k <- k + 1L
        }
      }
    }
  }

  class(out) <- c(switch(type,
                         ordered = "multi_ordered",
                         competing = "multi_competing"),
                  class(out))
  return(out)
}

#' Pull non-null variables from data
#'
#' Allows mutate syntax, e.g. factor(study) or sqrt(var)
#' Uses pull for variable names as strings or integer positions
#'
#' @param data data frame
#' @param var quosure (possibly NULL) for variable to pull
#'
#' @noRd
pull_non_null <- function(data, var) {
  var_null <- rlang::quo_is_missing(var) | rlang::quo_is_null(var)
  if (!var_null) {
    if (rlang::is_symbolic(rlang::quo_get_expr(var))) return(dplyr::pull(dplyr::transmute(data, {{ var }})))
    else return(dplyr::pull(data, {{ var }}))
  }
  else return(NULL)
}

#' Drop variables from original data set
#'
#' Allows for combination of pull and mutate semantics - doesn't drop mutated
#' results, handles integer column references
#'
#' @param data data frame after transformation
#' @param orig_data original data frame
#' @param var quosure for variable to drop
#'
#' @noRd
drop_original <- function(data, orig_data, var) {
  e_var <- rlang::quo_get_expr(var)
  if (rlang::is_symbol(e_var)) {
    # Character/bare column name, drop unless protected name (starts with dot)
    if (stringr::str_starts(rlang::as_name(e_var), "\\.")) return(data)
    else return(dplyr::select(data, - {{ var }}))
  } else if (rlang::is_integerish(e_var)) {
    # Integer column number, drop based on original location
    orig_var <- colnames(orig_data)[e_var]
    if (stringr::str_starts(orig_var, "\\.")) return(data)
    else return(dplyr::select(data, - {{ orig_var }}))
  } else {
    # Mutate expression, don't drop
    return(data)
  }
}

#' Get outcome type
#'
#' Determines outcome type based on which inputs are NA
#'
#' @noRd
get_outcome_type <- function(y, se, r, n, E, Surv) {
  o <- c()
  if (!is.null(y)) o <- c(o, "continuous")
  if (!is.null(r)) {
    if (inherits(r, "multi_ordered")) o <- c(o, "ordered")
    if (inherits(r, "multi_competing")) o <- c(o, "competing")
    if (!is.null(E)) o <- c(o, "rate")
    if (!is.null(n)) o <- c(o, "count")
    if (!inherits(r, c("multi_ordered", "multi_competing")) && is.null(n) && is.null(E)) o <- c(o, "binary")
  }
  if (!is.null(Surv)) o <- c(o, "survival")
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
#' @param append text to append to error message
#'
#' @noRd
check_outcome_continuous <- function(y, se = NULL, with_se = TRUE, append = NULL) {
  null_y <- is.null(y)
  null_se <- is.null(se)

  if (with_se) {
    if (!null_y && !null_se) {
      if (rlang::is_list(y) || !is.null(dim(y)))
        abort("Continuous outcome `y` must be a regular column (not a list or matrix column)")
      if (rlang::is_list(se) || !is.null(dim(se)))
        abort("Standard error `se` must be a regular column (not a list or matrix column)")
      if (!is.numeric(y)) abort(paste0("Continuous outcome `y` must be numeric", append))
      if (!is.numeric(se)) abort(paste0("Standard error `se` must be numeric", append))
      if (any(is.nan(se))) abort(paste0("Standard error `se` cannot be NaN", append))
      if (any(is.na(y))) abort(paste0("Continuous outcome `y` contains missing values", append))
      if (any(is.na(se))) abort(paste0("Standard error `se` contains missing values", append))
      if (any(is.infinite(se))) abort(paste0("Standard error `se` cannot be infinite", append))
      if (any(se <= 0)) abort(paste0("Standard errors must be positive", append))
    } else {
      if (!null_y) abort(paste0("Specify standard error `se` for continuous outcome `y`", append))
      if (!null_se) abort(paste0("Specify continuous outcome `y`", append))
    }
    invisible(list(y = y, se = se))
  } else {
    if (!null_y) {
      if (rlang::is_list(y) || !is.null(dim(y)))
        abort("Continuous outcome `y` must be a regular column (not a list or matrix column)")
      if (any(is.na(y))) abort(paste0("Continuous outcome `y` contains missing values", append))
      if (!is.numeric(y)) abort(paste0("Continuous outcome `y` must be numeric", append))
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
    if (rlang::is_list(n) || !is.null(dim(n)))
      abort("Denominator `n` must be a regular column (not a list or matrix column)")
    if (!is.numeric(n)) abort("Denominator `n` must be numeric")
    if (any(is.na(n))) abort("Denominator `n` contains missing values")
    if (any(n != trunc(n))) abort("Denominator `n` must be integer-valued")
    if (any(n <= 0)) abort("Denominator `n` must be greater than zero")
    if (null_r) abort("Specify outcome count `r`.")
  }

  if (!null_E) {
    if (rlang::is_list(E) || !is.null(dim(E)))
      abort("Time at risk `E` must be a regular column (not a list or matrix column)")
    if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
    if (any(is.na(E))) abort("Time at risk `E` contains missing values")
    if (any(E <= 0)) abort("Time at risk `E` must be positive")
    if (null_r) abort("Specify outcome count `r`.")
  }

  if (!null_r) {
    if (rlang::is_list(r) || !is.null(dim(r)))
      abort("Outcome count `r` must be a regular column (not a list or matrix column)")
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
      if (rlang::is_list(r) || !is.null(dim(r)))
        abort("Rate outcome count `r` must be a regular column (not a list or matrix column)")
      if (rlang::is_list(E) || !is.null(dim(E)))
        abort("Time at risk `E` must be a regular column (not a list or matrix column)")
      if (!is.numeric(E)) abort("Time at risk `E` must be numeric")
      if (any(is.na(E))) abort("Time at risk `E` contains missing values")
      if (any(E <= 0)) abort("Time at risk `E` must be positive")
      if (!is.numeric(r)) abort("Rate outcome count `r` must be numeric")
      if (any(is.na(r))) abort("Rate outcome count `r` contains missing values")
      if (any(r != trunc(r))) abort("Rate outcome count `r` must be non-negative integer")
      if (any(r < 0)) abort("Rate outcome count `r` must be non-negative integer")
    }
  } else if (!null_r) {
    if (rlang::is_list(r) || !is.null(dim(r)))
      abort("Binary outcome `r` must be a regular column (not a list or matrix column)")
    if (!is.numeric(r)) abort("Binary outcome `r` must be numeric")
    if (any(is.na(r))) abort("Binary outcome `r` contains missing values")
    if (any(! r %in% c(0, 1))) abort("Binary outcome `r` must equal 0 or 1")
  }

  invisible(list(r = r, E = E))
}

#' Check survival outcomes
#'
#' @param Surv vector
#'
#' @noRd
check_outcome_survival <- function(Surv) {
  if (!is.null(Surv)) {

    if (!survival::is.Surv(Surv)) abort("Survival outcome `Surv` must be a `Surv` object created using `Surv()`")

    stype <- attr(Surv, "type")
    allowed_stypes <- c("right", "left", "interval", "interval2", "counting")
    if (!stype %in% allowed_stypes)
      abort(glue::glue('Survival outcome `Surv` of type "{stype}" is not supported.\n',
                       "Supported types are ",
                       glue::glue_collapse(allowed_stypes, sep = ", ", last = " and "), "."))

    status <- Surv[, "status"]
    if (any(is.na(status))) abort("Survival outcome `Surv` contains missing event status values")
    if (!all(status %in% 0:3)) abort("Survival outcome `Surv` event status values must be 0, 1, 2, or 3")

    S <- get_Surv_data(Surv)
    if (any(is.na(S$time), is.na(S$start_time), is.na(S$delay_time))) abort("Survival outcome `Surv` contains missing times")
    if (any(is.infinite(S$time), is.infinite(S$start_time), is.infinite(S$delay_time))) abort("Survival outcome `Surv` contains infinite times")
    if (any(S$time <= 0)) abort("Survival outcome `Surv` must have strictly positive outcome times")
    if (any(S$start_time < 0, S$delay_time < 0)) abort("Survival outcome `Surv` must have non-negative start times")

  }

  invisible(list(Surv = Surv))
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
         ipd = c("continuous", NA)),
    list(agd_arm = c("ordered", NA),
         agd_contrast = c("continuous", NA),
         ipd = c("ordered", NA)),
    list(agd_arm = c("survival", NA),
         agd_contrast = c("continuous", NA),
         ipd = c("survival", NA))
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

#' Check sample size
#'
#' @param sample_size vector
#'
#' @noRd
check_sample_size <- function(sample_size) {
  if (rlang::is_list(sample_size) || !is.null(dim(sample_size)))
    abort("Sample size `sample_size` must be a regular column (not a list or matrix column)")
  if (!is.numeric(sample_size))
    abort("Sample size `sample_size` must be numeric")
  if (any(is.nan(sample_size)))
    abort("Sample size `sample_size` cannot be NaN")
  if (any(is.na(sample_size)))
    abort("Sample size `sample_size` contains missing values")
  if (any(sample_size != trunc(sample_size)))
    abort("Sample size `sample_size` must be integer-valued")
  if (any(sample_size <= 0))
    abort("Sample size `sample_size` must be greater than zero")
  if (any(is.infinite(sample_size)))
    abort("Sample size `sample_size` cannot be infinite")
}

#' Check treatment column
#'
#' @param trt Treatment vector
#'
#' @noRd
check_trt <- function(trt) {
  if (any(is.na(trt)))
    abort("`trt` cannot contain missing values")
  if (rlang::is_list(trt) || !is.null(dim(trt)))
    abort("`trt` must be a regular column (not a list or matrix column)")
}

#' Check study column
#'
#' @param study Treatment vector
#'
#' @noRd
check_study <- function(study) {
  if (any(is.na(study)))
    abort("`study` cannot contain missing values")
  if (rlang::is_list(study) || !is.null(dim(study)))
    abort("`study` must be a regular column (not a list or matrix column)")
}


#' Check treatment class coding
#'
#' @param trt_class Class vector
#' @param trt Treatment vector
#'
#' @noRd
check_trt_class <- function(trt_class, trt) {
  if (rlang::is_list(trt_class) || !is.null(dim(trt_class)))
    abort("`trt_class` must be a regular column (not a list or matrix column)")
  if (any(is.na(trt)))
    abort("`trt` cannot contain missing values")
  if (any(is.na(trt_class)))
    abort("`trt_class` cannot contain missing values")
  if (anyDuplicated(unique(cbind(trt, trt_class))[, "trt"]))
    abort("Treatment present in more than one class (check `trt` and `trt_class`)")
}

#' Checks for combining multinomial outcomes
#'
#' @param x List containing class `multi_*` outcome objects or the unclassed
#'   matrices within
#'
#' @noRd
check_multi_combine <- function(x) {
  # Remove any NULL elements
  x <- x[!purrr::map_lgl(x, is.null)]

  is_ordered <- purrr::map_lgl(x, ~inherits(., "multi_ordered"))
  is_competing <- purrr::map_lgl(x, ~inherits(., "multi_competing"))
  if (any(is_ordered) && any(is_competing))
    abort("Cannot combine ordered and competing multinomial outcomes.")

  x_u <- purrr::map(x, unclass)
  n_cat <- purrr::map_int(x_u, ncol)
  if (any(n_cat != n_cat[1]))
    abort("Cannot combine multinomial outcomes with different numbers of categories.")

  l_cat <- purrr::map(x_u, colnames)
  if (any(purrr::map_lgl(l_cat, ~any(. != l_cat[[1]]))))
    abort("Cannot combine multinomial outcomes with different category labels.")
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

#' Check whether AgD sample size columns are available
#'
#' @param network nma_data object
#'
#' @return logical TRUE/FALSE
#' @noRd
has_agd_sample_size <- function(network) {
  if (!inherits(network, "nma_data")) abort("Not nma_data object.")
  ss_a <- !has_agd_arm(network) || tibble::has_name(network$agd_arm, ".sample_size")
  ss_c <- !has_agd_contrast(network) || tibble::has_name(network$agd_contrast, ".sample_size")
  return(ss_a && ss_c)
}

#' Natural-order factors
#'
#' Produces factors with levels in natural sort order (i.e. 1 5 10 not 1 10 5)
#'
#' @noRd
nfactor <- function(x, ..., numeric = TRUE, resort = FALSE) {
  if (is.factor(x) && !resort) {
    return(x)
  } else {
    return(factor(x, levels = stringr::str_sort(unique(x), numeric = numeric), ...))
  }
}
