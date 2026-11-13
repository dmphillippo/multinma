#' @examples \donttest{
#' # Fitting a fixed effect model
#' smk_fit_FE <- nma(smk_net,
#'                   trt_effects = "fixed",
#'                   prior_intercept = normal(scale = 100),
#'                   prior_trt = normal(scale = 100))
#'
#' smk_fit_FE
#' }
#'
