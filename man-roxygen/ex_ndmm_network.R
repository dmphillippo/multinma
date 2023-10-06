#' @examples
#' # Set up newly-diagnosed multiple myeloma network
#'
#' head(ndmm_ipd)
#' head(ndmm_agd)
#'
#' ndmm_net <- combine_network(
#'   set_ipd(ndmm_ipd,
#'           study, trt,
#'           Surv = Surv(eventtime / 12, status)),
#'   set_agd_surv(ndmm_agd,
#'                study, trt,
#'                Surv = Surv(eventtime / 12, status),
#'                covariates = ndmm_agd_covs))
