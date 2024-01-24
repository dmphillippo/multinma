#' Working with 3D MCMC arrays
#'
#' 3D MCMC arrays (Iterations, Chains, Parameters) are produced by `as.array()`
#' methods applied to `stan_nma` or `nma_summary` objects.
#'
#' @rdname mcmc_array-class
#' @name mcmc_array-class
#' @aliases mcmc_array
#'
#' @examples
#' ## Smoking cessation
#' @template ex_smoking_nma_re_example
#' @examples \donttest{
#' # Working with arrays of posterior draws (as mcmc_array objects) is
#' # convenient when transforming parameters
#'
#' # Transforming log odds ratios to odds ratios
#' LOR_array <- as.array(relative_effects(smk_fit_RE))
#' OR_array <- exp(LOR_array)
#'
#' # mcmc_array objects can be summarised to produce a nma_summary object
#' smk_OR_RE <- summary(OR_array)
#'
#' # This can then be printed or plotted
#' smk_OR_RE
#' plot(smk_OR_RE, ref_line = 1)
#'
#' # Transforming heterogeneity SD to variance
#' tau_array <- as.array(smk_fit_RE, pars = "tau")
#' tausq_array <- tau_array^2
#'
#' # Correct parameter names
#' names(tausq_array) <- "tausq"
#'
#' # Summarise
#' summary(tausq_array)
#' }
NULL

#' @param x,object A 3D MCMC array of class `mcmc_array`
#' @param probs Numeric vector of quantiles of interest
#' @param ... Further arguments passed to other methods
#'
#' @rdname mcmc_array-class
#' @return The `summary()` method returns a [nma_summary] object, the `print()`
#'   method returns `x` invisibly. The `names()` method returns a character
#'   vector of parameter names, and `names()<-` returns the object with updated
#'   parameter names. The `plot()` method is a shortcut for
#'   `plot(summary(x), ...)`, passing all arguments on to [plot.nma_summary()].
#' @export
summary.mcmc_array <- function(object, ..., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  ss <- summary_mcmc_array(object, probs = probs)
  out <- list(summary = ss, sims = object)
  class(out) <- "nma_summary"
  return(out)
}

#' @rdname mcmc_array-class
#' @export
print.mcmc_array <- function(x, ...) {
  d <- dim(x)
  cglue("A MCMC array with {prod(d[1:2])} draws ({d[1]} iterations in {d[2]} chain{if (d[2] > 1) 's' else ''}) of {d[3]} parameter{if (d[3] > 1) 's' else ''}.")
  NextMethod(...)
  invisible(x)
}

#' @rdname mcmc_array-class
#' @export
plot.mcmc_array <- function(x, ...) {
  xsum <- list(summary = NULL, sims = x)
  class(xsum) <- "nma_summary"
  plot(xsum, ...)
}

#' @rdname mcmc_array-class
#' @export
names.mcmc_array <- function(x) {
  return(dimnames(x)[[3]])
}

#' @rdname mcmc_array-class
#' @param value Character vector of replacement parameter names
#' @export
`names<-.mcmc_array` <- function(x, value) {
  out <- x
  dimnames(out)[[3]] <- value
  return(out)
}
