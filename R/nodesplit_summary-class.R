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
#'   d_\mathrm{ind}}}
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
#' @param x A `nma_summary` object
#' @param ... Additional arguments passed on to other methods
#' @param digits Integer number of digits to display
#'
#' @rdname nodesplit_summary-methods
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
