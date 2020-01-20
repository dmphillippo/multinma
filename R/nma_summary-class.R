#' The `nma_summary` class
#'
#' The `nma_summary` class contains posterior summary statistics of model
#' parameters or other quantities of interest, and the draws used to obtain
#' these statistics.
#'
#' @rdname nma_summary-class
#' @name nma_summary-class
#' @aliases nma_summary
#'
#' @details Objects of class `nma_summary` have the following components:
#'   \describe{
#'   \item{summary}{A data frame containing the computed summary statistics}
#'   \item{sims}{A 3D array \[Iteration, Chain, Parameter\] of MCMC
#'   simulations}
#'   }
#'
NULL

#' Methods for `nma_summary` objects
#'
#' The `as.data.frame()`, `as_tibble()`, and `as.tibble()` methods return the
#' posterior summary statistics in a data frame or tibble. The `as.matrix()` and
#' `as.array()` methods return the posterior draws as a matrix or 3D array
#' (Iteration, Chain, Parameter).
#'
#'
#' @param x A `nma_summary` object
#' @param ... Additional arguments passed on to other methods
#'
#' @rdname nma_summary-methods
#'
#' @return
#' @export
#'
print.nma_summary <- function(x, ...) {
  # Format summaries nicely by study, if given
  print_study_block <- function(s, ...) {
    this_study <- unique(s$study)
    sec_header(this_study)
    s %>% dplyr::select(-.data$study) %>% print(...)
  }
  if (rlang::has_name(x, "study")) {
    by(x, x$study, print_study_block, ..., simplify = FALSE)
  } else {
    print(x$summary, ...)
  }
  invisible(x)
}

#' @rdname nma_summary-methods
#' @export
as.data.frame.nma_summary <- function(x, ...) {
  return(as.data.frame(x$summary, ...))
}

#' @rdname nma_summary-methods
#' @export
as.tibble.nma_summary <- function(x, ...) {
  return(x$summary)
}

#' @rdname nma_summary-methods
#' @export
as_tibble.nma_summary <- function(x, ...) {
  return(x$summary)
}

#' @rdname nma_summary-methods
#' @export
as.array.nma_summary <- function(x, ...) {
  return(x$sims)
}

#' @rdname nma_summary-methods
#' @export
as.matrix.nma_summary <- function(x, ...){
  # Follow approach in rstan:::as.matrix.stanfit
  a <- as.array(x)
  names_a <- dimnames(a)
  dim_a <- dim(a)
  dim(a) <- c(dim_a[1] * dim_a[2], dim_a[3])
  dimnames(a) <- names_a[-2]
  return(a)
}
