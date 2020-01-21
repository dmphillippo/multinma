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
#' @param digits Integer number of digits to display
#' @param pars Character vector of parameters to display in the printed summary
#' @param include Logical, are parameters named in `pars` included (`TRUE`) or excluded (`FALSE`)
#'
#' @rdname nma_summary-methods
#'
#' @return
#' @export
#'
print.nma_summary <- function(x, ..., digits = 2, pars, include) {
  # Set defaults for pars, include
  if (missing(include)) {
    include <- !missing(pars)
  } else {
    if (!is.logical(include) || length(include) > 1) abort("`include` should be TRUE or FALSE")
  }
  if (missing(pars)) {
    pars <- c("log_lik", "resdev", "fitted",
              "theta_bar_cum", "theta2_bar_cum",
              "mu", "delta", "lp__")
  } else {
    if (!is.character(pars)) abort("`pars` should be a character vector")
  }

  if (!is.numeric(digits) ||
      length(digits) > 1 ||
      trunc(digits) != digits) abort("`digits` must be a single integer")

  x_sum <- tibble::as_tibble(x) %>%
    dplyr::filter(
      if (include) grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
      else !grepl(paste0("^", pars, "(\\[.+\\])?$", collapse = "|"), .data$parameter)
      )

  # Get numeric columns for formatting, handling n_eff, Rhat separately
  num_col <- setdiff(names(x_sum)[purrr::map_lgl(x_sum, is.numeric)], c("n_eff", "Rhat"))

  x_sum <- dplyr::mutate_at(x_sum, num_col, ~round(., digits))

  if (rlang::has_name(x_sum, "n_eff")) x_sum$n_eff <- round(x_sum$n_eff, 0)
  if (rlang::has_name(x_sum, "Rhat")) x_sum$Rhat <- round(x_sum$Rhat, max(2, digits))

  x_sum <- tibble::column_to_rownames(x_sum, "parameter")

  # Format summaries nicely by study, if given
  print_study_block <- function(s, ...) {
    this_study <- unique(s$study)
    sec_header(this_study)
    s %>% dplyr::select(-.data$study) %>% print(...)
  }

  if (rlang::has_name(x_sum, "study")) {
    by(x_sum, x_sum$study, print_study_block, ..., simplify = FALSE)
  } else {
    print(x_sum, ...)
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
  return(tibble::as.tibble(x$summary, ...))
}

#' @rdname nma_summary-methods
#' @export
as_tibble.nma_summary <- function(x, ...) {
  return(tibble::as_tibble(x$summary, ...))
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
