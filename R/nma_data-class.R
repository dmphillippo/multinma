#' The nma_data class
#'
#' The `nma_data` class contains the data for a NMA in a standard format,
#' created using the functions [set_ipd()], [set_agd_arm()],
#' [set_agd_contrast()], or [combine_network()]. The sub-class `mlnmr_data` is
#' created by the function [add_integration()], and further contains numerical
#' integration points for the aggregate data.
#'
#' @rdname nma_data-class
#' @name nma_data-class
#' @aliases nma_data mlnmr_data mlnmr_data-class
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
#' Objects of class `mlnmr_data` additionally have components:
#'   \describe{
#'   \item{`n_int`}{number of numerical integration points}
#'   \item{`int_names`}{names of covariates with numerical integration points}
#'   }
#'
#' The `agd_arm` and `agd_contrast` tibbles have additional list columns with
#' prefix `.int_`, one for each covariate, which contain the numerical
#' integration points nested as length-`n_int` vectors within each row.
#'
NULL

#' Print `nma_data` objects
#'
#' @param x `nma_data` object
#' @param ... other options (not used)
#' @param n number of studies of each type to print
#'
#' @export
#'
#' @examples
print.nma_data <- function(x, ..., n = 10) {
  cwidth <- getOption("width")

  # Error if n is not an integer
  if (!is.numeric(n) | trunc(n) != n) abort("Argument `n` should be an integer")

  if (has_ipd(x)) {
    s_ipd <- x$ipd %>%
      dplyr::distinct(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise(Treatments = glue::glue("{dplyr::n()}: ",
                                              glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = .data$.study) %>%
      as.data.frame()
    n_ipd <- nrow(s_ipd)
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(x)) {
    s_agd_arm <- x$agd_arm %>%
      dplyr::distinct(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise(Treatments = glue::glue("{dplyr::n()}: ",
                                              glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = .data$.study) %>%
      as.data.frame()
    n_agd_arm <- nrow(s_agd_arm)
  } else {
    n_agd_arm <- 0
  }

  if (has_agd_contrast(x)) {
    s_agd_contrast <- x$agd_contrast %>%
      dplyr::distinct(.data$.study, .data$.trt, .data$.trt_b) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise(
        Treatments = glue::glue("{dplyr::n()}: ",
                     glue::glue_collapse(sort(unique(forcats::fct_c(.data$.trt, .data$.trt_b))),
                                         sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = .data$.study) %>%
      as.data.frame()
    n_agd_contrast <- nrow(s_agd_contrast)
  } else {
    n_agd_contrast <- 0
  }

  if (all(n_ipd == 0, n_agd_arm == 0, n_agd_contrast == 0)) {
    cglue("An empty `nma_data` object.")
  } else {
    cglue("A `nma_data` object with ", glue::glue_collapse(c(
      "{n_ipd} IPD stud{ifelse(n_ipd == 1, 'y', 'ies')}",
      "{n_agd_arm} AgD stud{ifelse(n_agd_arm == 1, 'y', 'ies')} (arm-based)",
      "{n_agd_contrast} AgD stud{ifelse(n_agd_contrast == 1, 'y', 'ies')} (contrast-based)"
    )[c(n_ipd > 0, n_agd_arm > 0, n_agd_contrast > 0)],
    last = ", and ", sep = ", "), ".")
  }
  cat("\n")

  if (n_ipd > 0) {
    sec_header("IPD studies")
    print(s_ipd[1:min(n_ipd, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_ipd > n) cglue(subtle(" ... plus {n_ipd - n} more studies"))
    cat("\n")
  }

  if (n_agd_arm > 0) {
    sec_header("AgD studies (arm-based)")
    print(s_agd_arm[1:min(n_agd_arm, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_agd_arm > n) cglue(subtle(" ... plus {n_agd_arm - n} more studies"))
    cat("\n")
  }

  if (n_agd_contrast > 0) {
    sec_header("AgD studies (contrast-based)")
    print(s_agd_contrast[1:min(n_agd_contrast, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_agd_contrast > n) cglue(subtle(" ... plus {n_agd_contrast - n} more studies"))
    cat("\n")
  }

  sec_header()
  cglue("Total number of treatments: {length(x$treatments)}")
  cglue("Total number of studies: {length(x$studies)}")
  cglue("Reference treatment is: {levels(x$treatments)[1]}")

  invisible(x)
}

#' @describeIn print.nma_data
#' @export
print.mlnmr_data <- function(x, ..., n = 10) {
  NextMethod()
  cat("\n")
  sec_header("Numerical integration")
  cglue("Numerical integration points available for {length(x$int_names)} covariates: ")
  cat(" ", x$int_names, "\n")
  cglue("Number of numerical integration points: {n_int}")
  invisible(x)
}

#' Make section headers for print.nma_class
#'
#' @param s section header string
#' @param width line width
#' @param sep separator between s and surrounding dashes
#'
#' @noRd
sec_header <- function(s = "",
                       width = min(80, getOption("width")),
                       sep = ifelse(nchar(s), " ", "")) {
  cat(subtle(strrep('-', width - nchar(s) - 2*nchar(sep))),
      crayon::bold(s),
      subtle("----"), "\n", sep = sep)
}

#' Cat out glue objects
#'
#' @param ... arguments passed to glue
#' @param sep string separator
#'
#' @noRd
cglue <- function(..., sep = "\n") {
  cat(glue::glue_data(.x = parent.frame(), ...), sep = sep)
}

#' Define crayon styles
#'
#' @noRd
subtle <- crayon::silver
emph_r <- crayon::red$bold
emph_g <- crayon::green$bold
