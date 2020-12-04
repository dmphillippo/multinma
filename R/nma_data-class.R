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
#'   \item{`classes`}{treatment class coding factor (same length as `treatments`
#'   for entire network)}
#'   \item{`studies`}{study coding factor for entire network}
#'   \item{`outcome`}{outcome type for each data source, named list}
#'   }
#'
#' The `agd_arm`, `agd_contrast`, and `ipd` components are
#' tibbles with the following columns:
#'   \describe{
#'   \item{`.study`}{study (as factor)}
#'   \item{`.trt`}{treatment (as factor)}
#'   \item{`.trtclass`}{treatment class (as factor), if specified}
#'   \item{`.y`}{continuous outcome}
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
#'   \item{`int_cor`}{correlation matrix for covariates used to generate
#'   numerical integration points}
#'   }
#'
#' The `agd_arm` and `agd_contrast` tibbles have additional list columns with
#' prefix `.int_`, one for each covariate, which contain the numerical
#' integration points nested as length-`n_int` vectors within each row.
#'
#' @template seealso_nma_data
NULL

#' Print `nma_data` objects
#'
#' Print details of networks stored as [nma_data] objects, as created by
#' [set_ipd()], [set_agd_arm()], [set_agd_contrast()], or [combine_network()].
#'
#' @param x `nma_data` object
#' @param ... other options (not used)
#' @param n number of studies of each type to print
#'
#' @export
print.nma_data <- function(x, ..., n = 10) {
  cwidth <- getOption("width")

  # Error if n is not an integer
  if (!rlang::is_scalar_integerish(n)) abort("Argument `n` should be an integer")

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
    cglue(" Outcome type: {x$outcome$ipd}",
          if (x$outcome$ipd == "ordered") " ({ncol(x$ipd$.r)} categories)" else "")
    # cat("\n")
  }

  if (n_agd_arm > 0) {
    sec_header("AgD studies (arm-based)")
    print(s_agd_arm[1:min(n_agd_arm, n), ], right = FALSE, row.names = FALSE, max = 9999L)
    if (n_agd_arm > n) cglue(subtle(" ... plus {n_agd_arm - n} more studies"))
    cat("\n")
    cglue(" Outcome type: {x$outcome$agd_arm}",
          if (x$outcome$agd_arm == "ordered") " ({ncol(x$agd_arm$.r)} categories)" else "")
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
  if (!is.null(x$classes)) {
    cglue("Total number of treatments: {length(x$treatments)}, in {nlevels(x$classes)} classes")
  } else {
    cglue("Total number of treatments: {length(x$treatments)}")
  }
  cglue("Total number of studies: {length(x$studies)}")
  cglue("Reference treatment is: {levels(x$treatments)[1]}")
  cglue("Network is {if (is_network_connected(x)) green('connected') else red('disconnected')}")

  invisible(x)
}

#' @rdname print.nma_data
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
  s <- as.character(s)
  cat(subtle(strrep('-', width - nchar(s) - 2*nchar(sep))),
      bold(s),
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
subtle <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::silver(...))
  else
    return(...)
}
bold <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::bold(...))
  else
    return(...)
}
red <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::red(...))
  else
    return(...)
}
green <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::green(...))
  else
    return(...)
}

#' Convert networks to graph objects
#'
#' The method `as.igraph()` converts `nma_data` objects into the form used by
#' the [igraph] package. The method `as_tbl_graph()` converts `nma_data` objects
#' into the form used by the [ggraph] and
#' \link[tidygraph:tidygraph-package]{tidygraph} packages.
#'
#' @param x An [nma_data] object to convert
#' @param ... Additional arguments
#' @param collapse Logical, collapse edges over studies? Default `TRUE`, only
#'   one edge is produced for each comparison (by IPD or AgD study type) with a
#'   `.nstudy` attribute giving the number of studies making that comparison. If
#'   `FALSE`, repeated edges are added for each study making the comparison.
#'
#' @return An `igraph` object for `as.igraph()`, a `tbl_graph` object for
#'   `as_tbl_graph()`.
#' @export
#'
#' @rdname graph_conversion
#'
#' @importFrom igraph as.igraph
#'
#' @template ex_smoking_network
#' @examples
#' # Convert to igraph object
#' igraph::as.igraph(smk_net)  # Edges combined by default
#' igraph::as.igraph(smk_net, collapse = FALSE)  # Without combining edges
#'
#' # Convert to tbl_graph object
#' tidygraph::as_tbl_graph(smk_net)  # Edges combined by default
#' tidygraph::as_tbl_graph(smk_net, collapse = FALSE)  # Without combining edges
as.igraph.nma_data <- function(x, ..., collapse = TRUE) {

  if (!rlang::is_bool(collapse))
    abort("`collapse` must be TRUE or FALSE.")

  if (has_ipd(x)) {
    e_ipd <- x$ipd %>%
      dplyr::distinct(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::group_modify(~make_contrasts(.x$.trt))

    if (collapse) {
      e_ipd <- e_ipd %>%
        dplyr::group_by(.data$.trt, .data$.trt_b) %>%
        dplyr::summarise(.nstudy = dplyr::n(), .type = "IPD")
    } else {
      e_ipd$.type <- "IPD"
    }

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
      dplyr::group_modify(~make_contrasts(.x$.trt))

    if (collapse) {
      e_agd <- e_agd %>%
        dplyr::group_by(.data$.trt, .data$.trt_b) %>%
        dplyr::summarise(.nstudy = dplyr::n(), .type = "AgD")
    } else {
      e_agd$.type <- "AgD"
    }

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
    dplyr::select(.data$from, .data$to, dplyr::everything())

  if (collapse) {
    e_all <- e_all %>%
      dplyr::mutate(.nstudy = dplyr::if_else(is.na(.data$.nstudy), 0L, .data$.nstudy))
  }

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

  if (!is.null(x$classes)) {
    v_all$.trtclass <- x$classes
  }

  g <- igraph::graph_from_data_frame(e_all, directed = FALSE, vertices = v_all)
  return(g)
}

#' @rdname graph_conversion
#'
#' @method as_tbl_graph nma_data
# Dynamically exported, see zzz.R
as_tbl_graph.nma_data <- function(x, ...) {
  require_pkg("tidygraph")
  return(tidygraph::as_tbl_graph(igraph::as.igraph(x, ...)))
}

#' Make treatment contrasts
#'
#' @param trt Vector of treatment codes
#'
#' @return
#' @noRd
make_contrasts <- function(trt) {
  contrs <- utils::combn(sort(trt), 2)
  return(dplyr::tibble(.trt = contrs[2,],
                       .trt_b = contrs[1,]))
}

#' Get default reference treatment
#'
#' @param network An `nma_data` object, as created by the functions `set_*()`,
#'   `combine_network()`, or `add_integration()`
#' @param ... Other arguments (unused)
#'
#' @return String
#' @noRd
get_default_trt_ref <- function(network, ...) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # g_c <- igraph::as.igraph(x)
  g_nc <- igraph::as.igraph(network, collapse = FALSE)

  nodes <-
    tibble::tibble(
      .trt = network$treatments,
      degree = igraph::degree(g_nc),
      betweenness = igraph::betweenness(g_nc)
    ) %>%
    dplyr::arrange(dplyr::desc(.data$degree),
                   dplyr::desc(.data$betweenness),
                   .data$.trt)

  return(as.character(nodes[[1, ".trt"]]))
}

#' Check network connectedness
#'
#' Check whether a network is connected - whether there is a path of study
#' evidence linking every pair of treatments in the network.
#'
#' @param network An `nma_data` object, as created by the functions `set_*()` or
#'   `combine_network()`.
#'
#' @return Logical `TRUE` or `FALSE`
#' @export
#'
#' @details Models will still run with disconnected networks. However, estimated
#'   relative effects between treatments across disconnected parts of the
#'   network will be entirely based on the prior distribution (typically very
#'   uncertain), as there is no information to update the prior distribution.
#'   Relative effects within each connected sub-network will be estimated as if
#'   each sub-network had been analysed separately.
#'
#' @examples ## Smoking cessation
#' @template ex_smoking_network
#' @examples is_network_connected(smk_net)  # TRUE, network is connected
#' @examples
#' ## A disconnected network
#' disc_net <- set_agd_arm(smoking[smoking$studyn %in% c(15, 21), ],
#'                         study = studyn,
#'                         trt = trtc,
#'                         r = r,
#'                         n = n)
#' is_network_connected(disc_net)  # FALSE, network is disconnected
#' disc_net
#' plot(disc_net)
is_network_connected <- function(network) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  return(igraph::is_connected(igraph::as.igraph(network)))
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
#' @param weight_nodes Weight nodes by the total sample size? Default is `FALSE`.
#' @param show_trt_class Colour treatment nodes by class, if `trt_class` is set?
#'   Default is `FALSE`.
#'
#' @details The default is equivalent to `layout = "linear"` and `circular =
#'   TRUE`, which places the treatment nodes on a circle in the order defined by
#'   the treatment factor variable. An alternative layout which may give good
#'   results for simple networks is `"sugiyama"`, which attempts to minimise the
#'   number of edge crossings.
#'
#'   `weight_nodes = TRUE` requires that sample sizes have been specified for
#'   any aggregate data in the network, using the `sample_size` option of
#'   `set_agd_*()`.
#'
#' @return A `ggplot` object, as produced by [ggraph()].
#' @export
#'
#' @examples ## Stroke prevention in atrial fibrillation
#' # Setting up the network
#' af_net <- set_agd_arm(atrial_fibrillation,
#'                       study = studyc,
#'                       trt = abbreviate(trtc, minlength = 3),
#'                       r = r,
#'                       n = n,
#'                       trt_class = trt_class)
#' af_net
#'
#' # Basic plot
#' plot(af_net)
#'
#' # Turn off weighting edges by number of studies
#' plot(af_net, weight_edges = FALSE)
#'
#' # Turn on weighting nodes by sample size
#' plot(af_net, weight_nodes = TRUE)
#'
#' # Colour treatment nodes by class
#' plot(af_net, weight_nodes = TRUE, show_trt_class = TRUE)
#'
#' # Output may be customised using standard ggplot commands
#' # For example, to display the legends below the plot:
#' plot(af_net, weight_nodes = TRUE, show_trt_class = TRUE) +
#'   ggplot2::theme(legend.position = "bottom",
#'                  legend.box = "vertical",
#'                  legend.margin = ggplot2::margin(0, 0, 0, 0),
#'                  legend.spacing = ggplot2::unit(0.5, "lines"))
#'
#' # Choosing a different ggraph layout, hiding some legends
#' plot(af_net, weight_nodes = TRUE, show_trt_class = TRUE,
#'      layout = "star") +
#'   ggplot2::guides(edge_width = "none", size = "none")
#'
plot.nma_data <- function(x, ..., layout, circular,
                          weight_edges = TRUE, weight_nodes = FALSE,
                          show_trt_class = FALSE) {
  if (missing(layout) && missing(circular)) {
    layout <- "linear"
    circular <- TRUE
  } else if (missing(layout)) {
    layout <- "linear"
  } else if (missing(circular)) {
    circular <- FALSE
  }

  if (!rlang::is_bool(weight_edges))
    abort("`weight_edges` must be TRUE or FALSE.")

  if (!rlang::is_bool(weight_nodes))
    abort("`weight_nodes` must be TRUE or FALSE.")

  if (weight_nodes && !has_agd_sample_size(x))
    abort(paste("AgD study sample sizes not specified in network, cannot weight nodes.",
                "Specify `sample_size` in set_agd_*(), or set weight_nodes = FALSE.", sep = "\n"))

  if (!rlang::is_bool(show_trt_class))
    abort("`show_trt_class` must be TRUE or FALSE.")

  if (show_trt_class && is.null(x$classes))
    abort(paste("Treatment classes not specified in network.",
                "Specify `trt_class` in set_*(), or set show_trt_class = FALSE.", sep = "\n"))

  dat_mixed <- has_ipd(x) && (has_agd_arm(x) || has_agd_contrast(x))
  g <- ggraph::ggraph(igraph::as.igraph(x), layout = layout, circular = circular, ...)

  if (weight_edges) {
    g <- g +
      ggraph::geom_edge_fan(ggplot2::aes(edge_width = .data$.nstudy,
                                       edge_colour = .data$.type),
                            lineend = "round") +
      ggraph::scale_edge_width_continuous("Number of studies",
                                          breaks = breaks_integer(positive = TRUE),
                                          limits = function(x) range(breaks_integer(positive = TRUE)(x)))
  } else {
    g <- g +
      ggraph::geom_edge_fan(ggplot2::aes(edge_colour = .data$.type),
                            edge_width = 1,
                            lineend = "round")
  }

  if (weight_nodes) {
    if (show_trt_class) {
      g <- g +
        ggraph::geom_node_point(ggplot2::aes(size = .data$.sample_size,
                                             fill = .data$.trtclass,
                                             colour = .data$.trtclass),
                                shape = 21)
    } else {
      g <- g +
        ggraph::geom_node_point(ggplot2::aes(size = .data$.sample_size),
                                fill = "grey90", colour = "grey60",
                                shape = 21)
    }
    g <- g +
      ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                             hjust = "outward", vjust = "outward") +
      ggplot2::scale_size_area("Total sample size", max_size = 12)
  } else {
    if (show_trt_class) {
      g <- g +
        ggraph::geom_node_label(ggplot2::aes(label = .data$name,
                                             fill = .data$.trtclass))
    } else {
      g <- g +
        ggraph::geom_node_label(ggplot2::aes(label = .data$name),
                                fill = "grey90")
    }
  }

  if (show_trt_class) {
    g <- g + ggplot2::scale_fill_discrete("Treatment Class", aesthetics = c("fill", "colour"))
  }

  g <- g +
    ggraph::scale_edge_colour_manual("Data", values = c(AgD = "#113259", IPD = "#55A480"),
                                     guide = if (dat_mixed) "legend" else FALSE) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::coord_fixed(clip = "off")
  return(g)
}

#' Automatic breaks for integer scales
#'
#' Compute automatic breaks for integer scales, making sure that the breaks are
#' also integers.
#'
#' @param prefer_n Set of preferred numbers of break points.
#' @param extend Extend the breaks beyond the range of the data? If `TRUE`
#'   ensures that breaks are equally spaced. Default `TRUE`.
#' @param positive Create breaks for strictly positive integers? Default `TRUE`.
#' @param ... Unused.
#'
#' @details If `extend = TRUE` (the default), breaks may lie outside of the
#'   range of the data, but are guaranteed to be equally spaced. If `extend =
#'   FALSE`, the limits of the breaks are equal to the smallest and largest
#'   values in the data, but within this range breaks may not be equally spaced.
#'
#' @return Returns a function to compute breaks, as required by the `breaks`
#'   argument to the `scale_*_continuous` functions in `ggplot2`.
#' @noRd
#'
#' @examples
#' breaks_integer()(1:9)
#' breaks_integer()(1:18) # extended, equally spaced
#' breaks_integer(extend = FALSE)(1:18) # un-extended, unequally spaced
#' breaks_integer(positive = FALSE)(-1:10) # allow negative values
breaks_integer <- function(prefer_n = c(5, 4, 3, 6), extend = TRUE, positive = TRUE, ...) {
  if (!rlang::is_integerish(prefer_n, finite = TRUE) || any(prefer_n < 2))
    abort("`prefer_n` must be a vector of integers greater than 1.")
  if (!rlang::is_bool(extend))
    abort("`extend` must be a logical value TRUE or FALSE.")
  if (!rlang::is_bool(positive))
    abort("`positive` must be a logical value TRUE or FALSE.")

  def_prefer_n <- prefer_n
  def_extend <- extend
  def_positive <- positive

  function(x, prefer_n = def_prefer_n, extend = def_extend, positive = def_positive, ...) {
    r <- range(x)
    l <- diff(r)

    if (l < max(prefer_n) - 1) {
      return(seq.int(r[1], r[2]))
    }

    n <- prefer_n[which.min(l %% (prefer_n - 1))]
    if (!extend || l %% (n - 1) == 0) {
      return(unique(round(seq.int(r[1], r[2], length.out = n))))
    } else {
      remainders <- l %% (prefer_n - 1)
      shifts <- prefer_n - 1 - remainders
      n <- prefer_n[which.min(shifts)]
      shift <- min(shifts)
      s_l <- if (positive) min(floor(shift / 2), r[1] - 1) else floor(shift / 2)
      s_u <- shift - s_l
      return(seq.int(r[1] - s_l, r[2] + s_u, length.out = n))
    }
  }
}
