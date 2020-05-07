#' @examples
#' # Fitting a fixed effect model
#' smk_fit_FE <- nma(smk_net,
#'               trt_effects = "fixed",
#'               prior_intercept = normal(scale = 100),
#'               prior_trt = normal(scale = 100))
#'
#' smk_fit_FE
#'
#' # Fitting a random effects model
#' smk_fit_RE <- nma(smk_net,
#'                   trt_effects = "random",
#'                   prior_intercept = normal(scale = 100),
#'                   prior_trt = normal(scale = 100),
#'                   prior_het = normal(scale = 5))
#'
#'  smk_fit_RE
#'
#' # Fitting an unrelated mean effects (inconsistency) model
#' smk_fit_RE_UME <- nma(smk_net,
#'                       consistency = "ume",
#'                       trt_effects = "random",
#'                       prior_intercept = normal(scale = 100),
#'                       prior_trt = normal(scale = 100),
#'                       prior_het = normal(scale = 5))
#'
#' smk_fit_RE_UME
#'
