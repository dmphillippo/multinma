#' @examples \donttest{
#' # Fitting a ML-NMR model.
#' # Specify a regression model to include effect modifier interactions for five
#' # covariates, along with main (prognostic) effects. We use a probit link and
#' # specify that the two-parameter Binomial approximation for the aggregate-level
#' # likelihood should be used. We set treatment-covariate interactions to be equal
#' # within each class. We narrow the possible range for random initial values with
#' # init_r = 0.1, since probit models in particular are often hard to initialise.
#' # Using the QR decomposition greatly improves sampling efficiency here, as is
#' # often the case for regression models.
#' pso_fit <- nma(pso_net, 
#'                trt_effects = "fixed",
#'                link = "probit",
#'                likelihood = "bernoulli2",
#'                regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
#'                class_interactions = "common",
#'                prior_intercept = normal(scale = 10),
#'                prior_trt = normal(scale = 10),
#'                prior_reg = normal(scale = 10),
#'                init_r = 0.1,
#'                QR = TRUE)
#' pso_fit
#' }
#'
