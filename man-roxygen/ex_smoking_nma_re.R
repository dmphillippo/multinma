#' @examples \donttest{
#' # Fitting a random effects model
#' smk_fit_RE <- nma(smk_net, 
#'                   trt_effects = "random",
#'                   prior_intercept = normal(scale = 100),
#'                   prior_trt = normal(scale = 100),
#'                   prior_het = normal(scale = 5))
#'
#' smk_fit_RE
#' }
#'
