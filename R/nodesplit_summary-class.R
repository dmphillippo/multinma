#' The `nodesplit_summary` class
#'
#' The `nodesplit_summary` class contains posterior summary statistics for
#' node-splitting models, as a result of calling `summary()` on a
#' `nma_nodesplit` or `nma_nodesplit_df` object.
#'
#' @rdname nodesplit_summary-class
#' @name nodesplit_summary-class
#' @aliases nodesplit_summary
#'
#' @details Objects of class `nodesplit_summary` are tibble data frames, with one row
#'   for each node-split comparison and columns:
#'   \describe{
#'   \item{`trt1`, `trt2`}{Treatments forming the comparison}
#'   \item{`summary`}{A list column containing [nma_summary] objects with the
#'   posterior summaries and draws for each of the node-splitting parameters}
#'   \item{`p_value`}{Bayesian p-value for inconsistency}
#'   \item{`dic`}{A list column containing [nma_dic] objects, giving the model
#'   fit statistics}
#'   }
#'
#'   The parameters included in `summary` are:
#'   \describe{
#'   \item{`d_net`}{Network estimate from the corresponding consistency model,
#'   if available}
#'   \item{`d_dir`}{Direct estimate from the node-splitting model}
#'   \item{`d_ind`}{Indirect estimate from the node-splitting model}
#'   \item{`omega`}{Inconsistency factor \eqn{\omega = d_\mathrm{dir} -
#'   d_\mathrm{ind}}{\omega = d_dir - d_ind}}
#'   \item{`tau`}{Heterogeneity standard deviation from the node-splitting
#'   model, if a random effects model was fitted}
#'   \item{`tau_consistency`}{Heterogeneity standard deviation from the
#'   corresponding consistency model, if available and if a random effects model
#'   was fitted}
#'   }
#'
#'
NULL

#' Methods for `nodesplit_summary` objects
#'
#' The `as.data.frame()`, `as_tibble()`, and `as.tibble()` methods return the
#' node-splitting summaries in a data frame or tibble.
#'
#' @param x A `nodesplit_summary` object
#' @param ... Additional arguments passed on to other methods
#' @param digits Integer number of digits to display
#'
#' @rdname nodesplit_summary-methods
#'
#' @return A `data.frame` for `as.data.frame()`, a `tbl_df` for `as.tibble()`
#'   and `as_tibble()`.
#'
#'   The `print()` method returns `x` invisibly.
#'
#' @seealso [plot.nodesplit_summary()]
#'
#' @export
#'
print.nodesplit_summary <- function(x, ..., digits = 2) {
  if (!rlang::is_scalar_integerish(digits)) abort("`digits` must be a single integer")

  n_ns <- nrow(x)

  cglue("Node-splitting model{if (n_ns > 1) 's' else ''} fitted for {n_ns} comparison{if (n_ns > 1) 's' else ''}.")

  for (i in 1:nrow(x)) {
    cglue("")
    if (n_ns > 1) {
      sec_header(glue::glue("Node-split {x$trt2[i]} vs. {x$trt1[i]}"))
      cglue("")
    }
    print(x$summary[[i]], ...)
    cglue("")
    print(x$dic[[i]])
    cglue("")
    cglue("Bayesian p-value: {format.pval(x$p_value[[i]], digits = digits, eps = 10^-digits)}")
  }

  invisible(x)
}

#' @export
#' @rdname nodesplit_summary-methods
#' @param nest Whether to return a nested tibble, with the full [nma_summary]
#'   and [nma_dic] objects, or to unnest their summaries, default `FALSE`
as_tibble.nodesplit_summary <- function(x, ..., nest = FALSE) {
  if (!rlang::is_bool(nest)) abort("`nest` must be a single logical value TRUE/FALSE.")

  if (nest) { # Return underlying nested tibble
    NextMethod(...)
  } else { # Unnest summary results
    out <- x
    out$summary <- purrr::map(out$summary, tibble::as_tibble)
    out <- tidyr::unnest(out, cols = "summary")
    out$dic <- purrr::map_dbl(out$dic, "dic")
    class(out) <- setdiff(class(out), "nodesplit_summary")
    return(out)
  }
}

#' @export
#' @rdname nodesplit_summary-methods
as.tibble.nodesplit_summary <- function(x, ..., nest = FALSE) {
  return(tibble::as_tibble(x, ..., nest = nest))
}

#' @export
#' @rdname nodesplit_summary-methods
as.data.frame.nodesplit_summary <- function(x, ...) {
  return(as.data.frame(tibble::as_tibble(x, ..., nest = FALSE)))
}
