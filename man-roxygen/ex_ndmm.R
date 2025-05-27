#' @examples \donttest{
#' # Fit Weibull (PH) model
#' ndmm_fit <- nma(ndmm_net, 
#'                 likelihood = "weibull",
#'                 prior_intercept = normal(scale = 100),
#'                 prior_trt = normal(scale = 10),
#'                 prior_aux = half_normal(scale = 10))
#'
#' ndmm_fit
#' }
