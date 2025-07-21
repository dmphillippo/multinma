#' MCMC Convergence Diagnostics for Stan NMA Models
#'
#' These functions provide convenient access to MCMC convergence diagnostics
#' for `stan_nma` objects, leveraging the `bayesplot` package for visualization
#' and Stan's built-in diagnostic capabilities.
#'
#' @name mcmc_diagnostics
#' @param x A `stan_nma` object created by [nma()]
#' @param pars Character vector of parameter names to include in diagnostics.
#'   If `NULL` (default), automatically selects a reasonable subset based on
#'   the model type. See Details for automatic parameter selection.
#' @param ... Additional arguments passed to the underlying `bayesplot` functions
#'
#' @details
#' **Automatic Parameter Selection:**
#'
#' When `pars = NULL`, the functions automatically select parameters based on
#' the model characteristics:
#' - Always included: treatment effects (`d`), intercept terms
#' - If random effects: heterogeneity parameters (`tau`)
#' - If regression: regression coefficients (`beta`)
#' - If auxiliary parameters: auxiliary effects (`aux`)
#' - Limits to maximum of 12 parameters for readability
#'
#' **Diagnostic Functions:**
#' - `mcmc_trace()`: Trace plots for visual assessment of chain mixing
#' - `mcmc_acf()`: Autocorrelation function plots to assess chain efficiency
#' - `mcmc_rhat()`: R-hat convergence diagnostics (should be < 1.1)
#' - `mcmc_neff()`: Effective sample size diagnostics
#' - `mcmc_diagnostics_plot()`: Combined diagnostic plot with multiple panels
#'
#' @return
#' - `mcmc_trace()`, `mcmc_acf()`: ggplot objects from bayesplot
#' - `mcmc_rhat()`, `mcmc_neff()`: ggplot objects showing diagnostic values
#' - `mcmc_diagnostics_plot()`: Combined ggplot with multiple diagnostic panels
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- nma(network, trt_effects = "random", ...)
#'
#' # Generate trace plots
#' mcmc_trace(fit)
#'
#' # Check R-hat diagnostics
#' mcmc_rhat(fit)
#'
#' # Comprehensive diagnostic plot
#' mcmc_diagnostics_plot(fit, type = c("trace", "rhat"))
#' }
#'
NULL

#' @rdname mcmc_diagnostics
#' @export
mcmc_trace <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("x must be a stan_nma object")
  }
  
  if (is.null(pars)) {
    pars <- get_default_parameters(x)
  }
  
  # Validate parameters exist
  pars <- validate_parameters(x, pars)
  
  # Extract posterior draws
  posterior_array <- rstan::extract(x$stanfit, pars = pars, permuted = FALSE)
  
  # Generate trace plot
  bayesplot::mcmc_trace(posterior_array, pars = pars, ...)
}

#' @rdname mcmc_diagnostics
#' @export
mcmc_acf <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("x must be a stan_nma object")
  }
  
  if (is.null(pars)) {
    pars <- get_default_parameters(x)
  }
  
  # Validate parameters exist
  pars <- validate_parameters(x, pars)
  
  # Extract posterior draws
  posterior_array <- rstan::extract(x$stanfit, pars = pars, permuted = FALSE)
  
  # Generate ACF plot
  bayesplot::mcmc_acf(posterior_array, pars = pars, ...)
}

#' @rdname mcmc_diagnostics
#' @export
mcmc_rhat <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("x must be a stan_nma object")
  }
  
  if (is.null(pars)) {
    pars <- get_default_parameters(x)
  }
  
  # Validate parameters exist
  pars <- validate_parameters(x, pars)
  
  # Extract R-hat values
  fit_summary <- rstan::summary(x$stanfit, pars = pars)$summary
  rhat_values <- fit_summary[, "Rhat"]
  
  # Generate R-hat plot
  bayesplot::mcmc_rhat(rhat_values, ...)
}

#' @rdname mcmc_diagnostics
#' @export
mcmc_neff <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("x must be a stan_nma object")
  }
  
  if (is.null(pars)) {
    pars <- get_default_parameters(x)
  }
  
  # Validate parameters exist
  pars <- validate_parameters(x, pars)
  
  # Extract effective sample size
  fit_summary <- rstan::summary(x$stanfit, pars = pars)$summary
  neff_values <- fit_summary[, "n_eff"]
  
  # Generate n_eff plot
  bayesplot::mcmc_neff(neff_values, ...)
}

#' @rdname mcmc_diagnostics
#' @param type Character vector specifying which diagnostic plots to include.
#'   Options: "trace", "acf", "rhat", "neff". Default includes all.
#' @export
mcmc_diagnostics_plot <- function(x, pars = NULL, 
                                  type = c("trace", "acf", "rhat", "neff"), ...) {
  if (!inherits(x, "stan_nma")) {
    abort("x must be a stan_nma object")
  }
  
  type <- match.arg(type, several.ok = TRUE)
  
  if (is.null(pars)) {
    pars <- get_default_parameters(x)
  }
  
  # Validate parameters exist
  pars <- validate_parameters(x, pars)
  
  # Generate individual plots
  plots <- list()
  
  if ("trace" %in% type) {
    plots$trace <- mcmc_trace(x, pars = pars, ...)
  }
  
  if ("acf" %in% type) {
    plots$acf <- mcmc_acf(x, pars = pars, ...)
  }
  
  if ("rhat" %in% type) {
    plots$rhat <- mcmc_rhat(x, pars = pars, ...)
  }
  
  if ("neff" %in% type) {
    plots$neff <- mcmc_neff(x, pars = pars, ...)
  }
  
  # Combine plots using patchwork
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(patchwork::wrap_plots(plots, ncol = 2))
  }
}

# Internal helper functions

#' Get default parameters for diagnostics based on model type
#' @param x stan_nma object
#' @return character vector of parameter names
#' @keywords internal
get_default_parameters <- function(x) {
  all_pars <- x$stanfit@model_pars
  
  # Start with core parameters
  default_pars <- character(0)
  
  # Always include treatment effects if present
  if ("d" %in% all_pars) {
    d_pars <- grep("^d\\[", rownames(rstan::summary(x$stanfit)$summary), value = TRUE)
    default_pars <- c(default_pars, head(d_pars, 6))  # Limit to first 6
  }
  
  # Include intercept/baseline parameters
  intercept_pars <- grep("^(alpha|mu|baseline)", all_pars, value = TRUE)
  default_pars <- c(default_pars, head(intercept_pars, 3))
  
  # Include heterogeneity parameters for random effects models
  if (x$trt_effects == "random") {
    het_pars <- grep("^(tau|sigma)", all_pars, value = TRUE)
    default_pars <- c(default_pars, head(het_pars, 2))
  }
  
  # Include regression coefficients if regression model
  if (!is.null(x$regression)) {
    beta_pars <- grep("^beta", all_pars, value = TRUE)
    default_pars <- c(default_pars, head(beta_pars, 3))
  }
  
  # Include auxiliary parameters if present
  aux_pars <- grep("^(aux|delta)", all_pars, value = TRUE)
  default_pars <- c(default_pars, head(aux_pars, 2))
  
  # Remove duplicates and limit total
  default_pars <- unique(default_pars)
  default_pars <- head(default_pars, 12)  # Maximum 12 parameters
  
  return(default_pars)
}

#' Validate that requested parameters exist in the model
#' @param x stan_nma object
#' @param pars character vector of parameter names
#' @return validated character vector of parameter names
#' @keywords internal
validate_parameters <- function(x, pars) {
  if (length(pars) == 0) {
    abort("No parameters specified or available for diagnostics")
  }
  
  # Get all available parameter names from summary
  all_available <- rownames(rstan::summary(x$stanfit)$summary)
  
  # Check which requested parameters exist
  missing_pars <- setdiff(pars, all_available)
  
  if (length(missing_pars) > 0) {
    warn(paste("The following parameters were not found in the model and will be ignored:",
               paste(missing_pars, collapse = ", ")))
    pars <- intersect(pars, all_available)
  }
  
  if (length(pars) == 0) {
    abort("None of the specified parameters were found in the model")
  }
  
  return(pars)
} 