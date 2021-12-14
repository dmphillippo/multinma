#' The nma_nodesplit class
#'
#' The `nma_nodesplit` and `nma_nodesplit_df` classes contains the results from
#' running a node-splitting model with the function [nma()].
#'
#' @rdname nma_nodesplit-class
#' @name nma_nodesplit-class
#' @aliases nma_nodesplit nma_nodesplit_df nma_nodesplit_df-class
#'
#' @details Objects of class `nma_nodesplit` inherit from the [stan_nma] class,
#'   and contain the results of fitting a single node-split model. They have one
#'   additional component, `nodesplit`, which gives the comparison that was
#'   node-split as a length 2 vector.
#'
#'   Objects of class `nma_nodesplit_df` are tibble data frames with one row
#'   for each node-split comparison and columns:
#'   \describe{
#'   \item{`trt1`, `trt2`}{Treatments forming the comparison}
#'   \item{`model`}{A list column containing the results of each model as a
#'   `nma_nodesplit` object}
#'   }
#'   Optionally, there will be an additional row for the consistency model if
#'   this was fitted (e.g. by `get_nodesplits(., include_consistency = TRUE)`)
#'   with `trt1` and `trt2` both `NA`.
#'
NULL

#' Print `nma_nodesplit_df` objects
#'
#' @param x A [nma_nodesplit_df] object
#' @param ... Further arguments passed to [print.stanfit()]
#'
#' @export
print.nma_nodesplit_df <- function(x, ...) {
  n_ns <- nrow(dplyr::filter(x, !is.na(.data$trt1) & !is.na(.data$trt2)))

  cglue("Node-splitting models fitted for {n_ns} comparisons.")
  cglue("To summarise these results, use `summary()`.")

  for (i in 1:nrow(x)) {
    cglue("")
    if (is.na(x$trt1[i]) && is.na(x$trt2[i])) { # Consistency model
      sec_header("Consistency model")
    } else {
      sec_header(glue::glue("Node-split {x$trt2[i]} vs. {x$trt1[i]}"))
    }
    cglue("")
    print(x$model[[i]], ...)
  }

  invisible(x)
}
