#' @examples \donttest{
#' # Fitting all possible node-splitting models
#' smk_fit_RE_nodesplit <- nma(smk_net, \dontshow{refresh = if (interactive()) 200 else 0,}
#'                             consistency = "nodesplit",
#'                             trt_effects = "random",
#'                             prior_intercept = normal(scale = 100),
#'                             prior_trt = normal(scale = 100),
#'                             prior_het = normal(scale = 5))
#' }
#'
