#' The nma_data class
#'
#' The `nma_data` class contains the data for a NMA in a standard format,
#' created using the functions [set_ipd()], [set_agd_arm()],
#' [set_agd_contrast()], [set_agd_surv()], or [combine_network()]. The sub-class
#' `mlnmr_data` is created by the function [add_integration()], and further
#' contains numerical integration points for the aggregate data.
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
#'   \item{`.Surv`}{survival outcome of type [`Surv`] (time-to-event), nested by
#'   study arm}
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
#' [set_ipd()], [set_agd_arm()], [set_agd_contrast()], [set_agd_surv()], or
#' [combine_network()].
#'
#' @param x `nma_data` object
#' @param ... other options (not used)
#' @param n number of studies of each type to print
#'
#' @return `x` is returned invisibly.
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
      dplyr::summarise("Treatment arms" = glue::glue("{dplyr::n()}: ",
                                              glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = ".study") %>%
      as.data.frame()
    n_ipd <- nrow(s_ipd)
  } else {
    n_ipd <- 0
  }

  if (has_agd_arm(x)) {
    s_agd_arm <- x$agd_arm %>%
      dplyr::arrange(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise("Treatment arms" = glue::glue("{dplyr::n()}: ",
                                              glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = ".study") %>%
      as.data.frame()
    n_agd_arm <- nrow(s_agd_arm)
  } else {
    n_agd_arm <- 0
  }

  if (has_agd_contrast(x)) {
    s_agd_contrast <- x$agd_contrast %>%
      dplyr::arrange(.data$.study, .data$.trt) %>%
      dplyr::group_by(.data$.study) %>%
      dplyr::summarise("Treatment arms" = glue::glue("{dplyr::n()}: ",
                                               glue::glue_collapse(.data$.trt, sep = " | ", width = 0.8*cwidth))) %>%
      dplyr::rename(Study = ".study") %>%
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
    return(paste0(...))
}
bold <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::bold(...))
  else
    return(paste0(...))
}
red <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::red(...))
  else
    return(paste0(...))
}
green <- function(...) {
  if (require_pkg("crayon", error = FALSE))
    return(crayon::green(...))
  else
    return(paste0(...))
}

#' Convert networks to graph objects
#'
#' The method `as.igraph()` converts `nma_data` objects into the form used by
#' the \link[igraph:igraph-package]{igraph} package. The method `as_tbl_graph()`
#' converts `nma_data` objects into the form used by the
#' \link[ggraph:ggraph-package]{ggraph} and
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
        dplyr::summarise(.nstudy = dplyr::n_distinct(.data$.study), .type = "IPD")
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
        dplyr::summarise(.nstudy = dplyr::n_distinct(.data$.study), .type = "AgD")
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
    dplyr::rename(from = ".trt_b", to = ".trt") %>%
    dplyr::select("from", "to", dplyr::everything())

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

  if (!is.null(x$classes)) v_all$.trtclass <- x$classes

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
#' @noRd
make_contrasts <- function(trt) {
  if (length(trt) < 2) {
    return(dplyr::tibble(.trt = trt, .trt_b = trt))
  }

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

  # Prefer longer follow-up for survival outcomes
  surv_dat <- dplyr::bind_rows(
    if (identical(network$outcome$ipd, "survival")) dplyr::select(network$ipd, ".trt", ".Surv") else NULL,
    if (identical(network$outcome$agd_arm, "survival")) tidyr::unnest(network$agd_arm, cols = ".Surv") %>% dplyr::select(".trt", ".Surv") else NULL
  )

  if (nrow(surv_dat) < 1) {
    max_fu <- 0
  } else {
    max_fu <- dplyr::mutate(surv_dat, !!! get_Surv_data(surv_dat$.Surv)) %>%
      dplyr::group_by(.data$.trt) %>%
      dplyr::summarise(max_fu = max(.data$time)) %>%
      dplyr::pull("max_fu")
  }

  g_nc <- igraph::as.igraph(network, collapse = FALSE)

  nodes <-
    tibble::tibble(
      .trt = network$treatments,
      degree = igraph::degree(g_nc),
      betweenness = igraph::betweenness(g_nc),
      max_fu = max_fu
    ) %>%
    dplyr::arrange(dplyr::desc(.data$degree),
                   dplyr::desc(.data$betweenness),
                   dplyr::desc(.data$max_fu),
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

#' Direct and indirect evidence
#'
#' Determine whether two treatments in a network are connected by direct and/or
#' indirect evidence, and generate a list of comparisons with both direct and
#' indirect evidence (i.e. potential inconsistency) for node-splitting.
#'
#' @param network An `nma_data` object, as created by the functions `set_*()` or
#'   `combine_network()`.
#' @param include_consistency Logical, whether to include a row of `NA`s to
#'   indicate that a consistency model (i.e. a model with no node-splitting)
#'   should also be fitted by the [nma()] function. Default is `FALSE` when
#'   calling `get_nodesplits()` by hand, and [nma()] sets this to `TRUE` by
#'   default.
#'
#' @details The list of comparisons for node-splitting is generated following
#'   the algorithm of \insertCite{Valkenhoef2016;textual}{multinma}. A
#'   comparison between two treatments has the potential for inconsistency, and
#'   is thus considered for node-splitting, if the comparison has both direct
#'   evidence and *independent* indirect evidence.
#'
#'   The notion of independent indirect evidence is necessary when multi-arm
#'   trials are present, since by design these trials are internally consistent.
#'   A comparison between two treatments has independent indirect evidence if,
#'   after removing all studies comparing the two treatments from the network,
#'   the two treatments are still connected by a path of evidence. This is the
#'   criterion considered by the `has_indirect()` function.
#'
#' @return For `has_direct()` and `has_indirect()`, a single logical value. For
#'   `get_nodesplits()`, a data frame with two columns giving the comparisons
#'   for node-splitting.
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' # Parkinsons example
#' park_net <- set_agd_arm(parkinsons,
#'                         study = studyn,
#'                         trt = trtn,
#'                         y = y,
#'                         se = se,
#'                         trt_ref = 1)
#'
#' # View the network plot
#' plot(park_net)
#'
#' # The 4 vs. 5 comparison is a spur on the network
#' has_direct(park_net, 4, 5)
#' has_indirect(park_net, 4, 5)
#'
#' # 1 and 5 are not directly connected
#' has_direct(park_net, 1, 5)
#' has_indirect(park_net, 1, 5)
#'
#' # The 1 vs. 2 comparison does not have independent indirect evidence, since
#' # the 1-2-4 loop is a multi-arm study
#' has_indirect(park_net, 1, 2)
#'
#' # Get a list of comparisons with potential inconsistency for node-splitting
#' get_nodesplits(park_net)
#'
#' # See van Valkenhoef (2016) for a discussion of this example

get_nodesplits <- function(network, include_consistency = FALSE) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # Check other arguments
  if (!rlang::is_bool(include_consistency))
    abort("`include_consistency` must be a logical value (TRUE or FALSE).")

  # Determine which contrasts to split using algorithm of van Valkenhoef,
  # i.e. having both direct and *independent* indirect evidence
  comparisons <- igraph::as_edgelist(igraph::as.igraph(network))
  colnames(comparisons) <- c("trt1", "trt2")

  out <- dplyr::as_tibble(comparisons) %>%
    # Remove edges of treatment against itself (from studies with multiple arms of the same treatment)
    dplyr::filter(.data$trt1 != .data$trt2) %>%
    dplyr::mutate(trt1 = factor(.data$trt1, levels = levels(network$treatments)),
                  trt2 = factor(.data$trt2, levels = levels(network$treatments))) %>%
    dplyr::arrange(.data$trt1, .data$trt2) %>%
    dplyr::rowwise() %>%
    dplyr::filter(has_direct(network, .data$trt1, .data$trt2) &&
                    has_indirect(network, .data$trt1, .data$trt2)) %>%
    dplyr::ungroup()

  # Add an NA row for the consistency model if include_consistency = TRUE
  if (include_consistency) {
    out <- dplyr::add_row(out, trt1 = NA, trt2 = NA, .before = 1)
  }

  return(out)
}

#' @param trt1,trt2 Treatments, each as a single integer, string, or factor
#' @export
#' @rdname get_nodesplits
has_direct <- function(network, trt1, trt2) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # Check treatments
  if (!rlang::is_scalar_atomic(trt1))
    abort("`trt1` should be a single integer, string, or factor, naming a treatment.")
  if (!rlang::is_scalar_atomic(trt2))
    abort("`trt2` should be a single integer, string, or factor, naming a treatment.")

  trt1 <- as.character(trt1)
  trt2 <- as.character(trt2)
  lvls_trt <- levels(network$treatments)

  if (! trt1 %in% lvls_trt)
    abort(sprintf("`trt1` does not match a treatment in the network.\nSuitable values are: %s",
                  ifelse(length(lvls_trt) <= 5,
                         paste0(lvls_trt, collapse = ", "),
                         paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

  if (! trt2 %in% lvls_trt)
    abort(sprintf("`trt2` does not match a treatment in the network.\nSuitable values are: %s",
                  ifelse(length(lvls_trt) <= 5,
                         paste0(lvls_trt, collapse = ", "),
                         paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

  if (trt1 == trt2)
    abort("`trt1` and `trt2` cannot be the same treatment.")

  # Convert to igraph and return whether adjacent nodes or not
  return(igraph::are_adjacent(igraph::as.igraph(network), trt1, trt2))
}

#' @rdname get_nodesplits
#' @export
has_indirect <- function(network, trt1, trt2) {

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("`network` must be an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  # Check treatments
  if (!rlang::is_scalar_atomic(trt1))
    abort("`trt1` should be a single integer, string, or factor, naming a treatment.")
  if (!rlang::is_scalar_atomic(trt2))
    abort("`trt2` should be a single integer, string, or factor, naming a treatment.")

  trt1 <- as.character(trt1)
  trt2 <- as.character(trt2)
  lvls_trt <- levels(network$treatments)

  if (! trt1 %in% lvls_trt)
    abort(sprintf("`trt1` does not match a treatment in the network.\nSuitable values are: %s",
                  ifelse(length(lvls_trt) <= 5,
                         paste0(lvls_trt, collapse = ", "),
                         paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

  if (! trt2 %in% lvls_trt)
    abort(sprintf("`trt2` does not match a treatment in the network.\nSuitable values are: %s",
                  ifelse(length(lvls_trt) <= 5,
                         paste0(lvls_trt, collapse = ", "),
                         paste0(paste0(lvls_trt[1:5], collapse = ", "), ", ..."))))

  if (trt1 == trt2)
    abort("`trt1` and `trt2` cannot be the same treatment.")


  # Create network with studies on both trt1 and trt2 removed
  studies <- dplyr::bind_rows(
    if (identical(network$outcome$agd_arm, "survival")) dplyr::select(network$agd_arm, -.data$.Surv) else network$agd_arm,
    network$agd_contrast,
    if (identical(network$outcome$ipd, "survival")) dplyr::select(network$ipd, -.data$.Surv) else network$ipd
  ) %>%
    dplyr::distinct(.data$.study, .data$.trt)

  ind_studies <- studies %>%
    dplyr::group_by(.data$.study) %>%
    dplyr::filter(! all(c(trt1, trt2) %in% .data$.trt))

  # Catch case where the reduced network is empty
  if (nrow(ind_studies) == 0) return(FALSE)

  ind_e <- dplyr::group_modify(ind_studies, ~make_contrasts(.x$.trt)) %>%
    dplyr::transmute(from = .data$.trt_b,
                     to = .data$.trt,
                     .data$.study)
  ind_v <- dplyr::ungroup(studies) %>% dplyr::distinct(.data$.trt)

  g <- igraph::graph_from_data_frame(ind_e, directed = FALSE, vertices = ind_v)

  # There is indirect evidence if the treatments are still connected in this
  # reduced network
  return(!is.infinite(igraph::distances(g, trt1, trt2)[1, 1]))
}

#' Network plots
#'
#' Create a network plot from a `nma_data` network object.
#'
#' @param x A [nma_data] object to plot
#' @param ... Additional arguments passed to
#'   \code{\link[ggraph:ggraph]{ggraph()}} and on to the layout function
#' @param layout The type of layout to create. Any layout accepted by
#'   \code{\link[ggraph:ggraph]{ggraph()}} may be used, including all of the
#'   layout functions provided by \link[igraph:igraph-package]{igraph}.
#' @param circular Whether to use a circular representation. See
#'   \code{\link[ggraph:ggraph]{ggraph()}}.
#' @param weight_edges Weight edges by the number of studies? Default is `TRUE`.
#' @param weight_nodes Weight nodes by the total sample size? Default is
#'   `FALSE`.
#' @param show_trt_class Colour treatment nodes by class, if `trt_class` is set?
#'   Default is `FALSE`.
#' @param level Display network at the `"treatment"` (default) or `"class"` level.
#' @param nudge Numeric value to nudge the treatment labels away from the nodes
#'   when `weight_nodes = TRUE`. Default is `0` (no adjustment to label
#'   position). A small value like `0.1` is usually sufficient.
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
#' @return A `ggplot` object, as produced by
#'   \code{\link[ggraph:ggraph]{ggraph()}}.
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
#' # Nudge the treatment labels away from the nodes
#' plot(af_net, weight_nodes = TRUE, show_trt_class = TRUE, nudge = 0.1)
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
                          weight_edges = TRUE,
                          weight_nodes = FALSE,
                          show_trt_class = FALSE,
                          level = c("treatment", "class"),
                          nudge = 0) {
  level <- rlang::arg_match(level)
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

  if(level == "class" && is.null(x$classes))
    abort(paste("Treatment classes not specified in network.",
                'Specify `trt_class` in set_*(), or use level = "treatment"', sep = "\n"))

  if (!rlang::is_double(nudge, n = 1, finite = TRUE))
    abort("`nudge` must be a single numeric value")

  # network at the class level
  if (level == "class"){
    if(has_agd_contrast(x)){
      x$agd_contrast$.trt <- x$agd_contrast$.trtclass
    }
    if(has_agd_arm(x)){
      x$agd_arm$.trt <- x$agd_arm$.trtclass
    }
    if(has_ipd(x)){
      x$ipd$.trt <- x$ipd$.trtclass
    }
    x$classes <- forcats::fct_unique(x$classes)
    x$treatments <- x$classes
  }

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

    # Calculate nudge positions
    if (nudge == 0) {
      pos <- ggplot2::position_identity()
    } else {
      if (circular || (rlang::is_string(layout) && layout %in% c("circle", "star"))) {
        pos <- ggplot2::position_nudge(x = nudge * g$data$x,
                                       y = nudge * g$data$y)
      } else {
        pos <- ggplot2::position_nudge(x = nudge * sign(g$data$x - mean(g$data$x)),
                                       y = nudge * sign(g$data$y - mean(g$data$y)))
      }
    }

    g <- g +
      ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                             hjust = "outward", vjust = "outward",
                             position = pos) +
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
                                     guide = if (dat_mixed) "legend" else "none") +
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
