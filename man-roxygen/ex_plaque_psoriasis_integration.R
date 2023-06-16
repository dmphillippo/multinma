#' @examples
#' # Add integration points to the network
#' pso_net <- add_integration(pso_net,
#'   durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
#'   prevsys = distr(qbern, prob = prevsys),
#'   bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
#'   weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
#'   psa = distr(qbern, prob = psa),
#'   n_int = 64)
#'
