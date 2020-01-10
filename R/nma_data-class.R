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
#'   \item{`outcome`}{outcome type for each data source, named list}
#'   }
#'
#' The `agd_arm`, `agd_contrast`, and `ipd` components are
#' tibbles with the following columns:
#'   \describe{
#'   \item{`.study`}{study (as factor)}
#'   \item{`.trt`}{treatment (as factor)}
#'   \item{`.y`}{outcome (continuous)}
#'   \item{`.se`}{standard error (continuous)}
#'   \item{`.r`}{event count (discrete)}
#'   \item{`.n`}{event count denominator (discrete, `agd_arm` only)}
#'   \item{`.E`}{time at risk (discrete)}
##'   \item{`.surv`}{event/censoring time, of type `Surv` (time-to-event)}
#'   \item{`.sample_size`}{sample size (`agd_*` only)}
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
      dplyr::distinct(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise(Treatments = glue::glue("{dplyr::n()}: ",
                                               glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = .data$.study) %>%
      as.data.frame()
    n_agd_contrast <- nrow(s_agd_contrast)
  } else {
    n_agd_contrast <- 0
  }

  if (all(n_ipd == 0, n_agd_arm == 0, n_agd_contrast == 0)) {
    cglue("An empty network.")
  } else {
    cglue("A network with ", glue::glue_collapse(c(
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
    cglue(" Outcome type: {x$outcome$ipd}")
    # cat("\n")
  }

  if (n_agd_arm > 0) {
    sec_header("AgD studies (arm-based)")
    print(s_agd_arm[1:min(n_agd_arm, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_agd_arm > n) cglue(subtle(" ... plus {n_agd_arm - n} more studies"))
    cat("\n")
    cglue(" Outcome type: {x$outcome$agd_arm}")
    # cat("\n")
  }

  if (n_agd_contrast > 0) {
    sec_header("AgD studies (contrast-based)")
    print(s_agd_contrast[1:min(n_agd_contrast, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_agd_contrast > n) cglue(subtle(" ... plus {n_agd_contrast - n} more studies"))
    cat("\n")
    cglue(" Outcome type: {x$outcome$agd_contrast}")
    # cat("\n")
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
  cat(x$int_names, fill = TRUE, labels = " ")
  cglue("Number of numerical integration points: {x$n_int}")
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
subtle <- function(...) crayon::silver(...)
emph_r <- function(...) crayon::red$bold(...)
emph_g <- function(...) crayon::green$bold(...)

#' Convert networks to graph objects
#'
#' The method `as.igraph()` converts `nma_data` objects into the form used by
#' the [igraph] package. The method `as_tbl_graph()` converts `nma_data` objects
#' into the form used by the [ggraph] and [tidygraph] packages.
#'
#' @param x An [nma_data] object to convert
#' @param ... Additional arguments
#'
#' @return
#' @export
#'
#' @rdname graph_conversion
#'
#' @importFrom igraph as.igraph
#'
#' @examples
as.igraph.nma_data <- function(x, ...) {

  if (has_ipd(x)) {
    e_ipd <- x$ipd %>%
      dplyr::distinct(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::group_modify(~make_contrasts(.x$.trt)) %>%
      dplyr::group_by(.data$.trt, .data$.trt_b) %>%
      dplyr::summarise(.nstudy = dplyr::n(), .type = "IPD")

    v_ipd <- x$ipd %>%
      dplyr::group_by(.data$.trt) %>%
      dplyr::summarise(.sample_size = dplyr::n())
  } else {
    e_ipd <- v_ipd <- tibble::tibble()
  }

  if (has_agd_arm(x) || has_agd_contrast(x)) {
    agd_all <- dplyr::bind_rows(x$agd_arm, x$agd_contrast)
    e_agd <- agd_all %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::group_modify(~make_contrasts(.x$.trt)) %>%
      dplyr::group_by(.data$.trt, .data$.trt_b) %>%
      dplyr::summarise(.nstudy = dplyr::n(), .type = "AgD")

    if (has_agd_sample_size(x)) {
      v_agd <- agd_all %>%
        dplyr::group_by(.data$.trt) %>%
        dplyr::summarise(.sample_size = sum(.data$.sample_size))
    } else {
      v_agd <- agd_all %>% dplyr::distinct(.data$.trt)
    }
  } else {
    e_agd <- v_agd <- tibble::tibble()
  }

  e_all <- dplyr::bind_rows(e_ipd, e_agd) %>%
    dplyr::rename(from = .data$.trt_b, to = .data$.trt) %>%
    dplyr::mutate(.nstudy = dplyr::if_else(is.na(.data$.nstudy), 0L, .data$.nstudy))

  if (has_agd_sample_size(x)) {
    v_all <- dplyr::bind_rows(v_ipd, v_agd) %>%
      dplyr::group_by(.data$.trt) %>%
      dplyr::summarise(.sample_size = sum(.data$.sample_size, na.rm = TRUE)) %>%
      dplyr::arrange(.data$.trt)
  } else {
    v_all <- dplyr::bind_rows(v_ipd, v_agd) %>%
      dplyr::distinct(.data$.trt) %>%
      dplyr::arrange(.data$.trt)
  }

  g <- igraph::graph_from_data_frame(e_all, directed = FALSE, vertices = v_all)
  return(g)
}

#' @return
#' @export
#'
#' @rdname graph_conversion
#'
#' @importFrom tidygraph as_tbl_graph
#'
#' @examples
as_tbl_graph.nma_data <- function(x, ...) {
  return(tidygraph::as_tbl_graph(igraph::as.igraph(x, ...)))
}

#' Make treatment contrasts
#'
#' @param trt Vector of treatment codes
#'
#' @return
#' @noRd
make_contrasts <- function(trt) {
  contrs <- utils::combn(trt, 2)
  return(dplyr::tibble(.trt = contrs[2,],
                       .trt_b = contrs[1,]))
}


#' Network plots
#'
#' Create a network plot from a `nma_data` network object.
#'
#' @param x A [nma_data] object to plot
#' @param ... Additional arguments passed to [ggraph()] and on to the layout
#'   function
#' @param layout The type of layout to create. Any layout accepted by [ggraph()]
#'   may be used, including all of the layout functions provided by [igraph].
#' @param circular Whether to use a circular representation. See [ggraph()].
#' @param weight_edges Weight edges by the number of studies? Default is `TRUE`.
#'
#' @details The default is equivalent to `layout = "linear"` and `circular =
#'   TRUE`, which places the treatment nodes on a circle in the order defined by
#'   the treatment factor variable. An alternative layout which may give good
#'   results for simple networks is `"sugiyama"`, which attempts to minimise the
#'   number of edge crossings.
#'
#' @return
#' @export
#'
#' @examples
plot.nma_data <- function(x, ..., layout, circular, weight_edges = TRUE) {
  if (missing(layout) && missing(circular)) {
    layout <- "linear"
    circular <- TRUE
  } else if (missing(layout)) {
    layout <- "linear"
  } else if (missing(circular)) {
    circular <- FALSE
  }

  if (!is.logical(weight_edges) || length(weight_edges) > 1)
    abort("`weight_edges` must be TRUE or FALSE.")

  dat_mixed <- has_ipd(x) && (has_agd_arm(x) || has_agd_contrast(x))
  g <-
    ggraph::ggraph(x, layout = layout, circular = circular, ...) +
    {if (weight_edges) {
      ggraph::geom_edge_fan(ggplot2::aes(edge_width = .data$.nstudy,
                                         edge_colour = .data$.type),
                            lineend = "round")
    } else {
      ggraph::geom_edge_fan(ggplot2::aes(edge_colour = .data$.type),
                            edge_width = 1,
                            lineend = "round")
    }} +
    ggraph::geom_node_label(ggplot2::aes(label = .data$name), fill = "grey90") +
    ggraph::scale_edge_colour_manual("Data", values = c(AgD = "#113259", IPD = "#55A480"),
                                     guide = if (dat_mixed) "legend" else FALSE) +
    {if (weight_edges) ggraph::scale_edge_width_continuous("Number of studies")} +
    ggraph::theme_graph(base_family = "") +
    ggplot2::coord_cartesian(clip = "off")
  return(g)
}
