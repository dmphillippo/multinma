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
