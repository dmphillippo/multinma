#' Example smoking FE NMA
#'
#' Smoking FE NMA for use in examples.
#'
#' @name example_smk_fe
#' @rdname aa_example_smk_fe
#' @description Calling `example("example_smk_fe")` will run a fixed effects
#'   NMA model with the smoking cessation data, using the code in the Examples
#'   section below. The resulting `stan_nma` object `smk_fit_FE` will then be
#'   available in the global environment.
#'
#' @keywords examples
#' @template ex_smoking_network
#' @template ex_smoking_nma_fe
#' @examples \dontshow{
#' if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) {
#'   assign("smk_net", smk_net, .GlobalEnv)
#'   assign("smk_fit_FE", smk_fit_FE, .GlobalEnv)
#' }
#' }
NULL

#' Example smoking RE NMA
#'
#' Smoking RE NMA for use in examples.
#'
#' @name example_smk_re
#' @rdname aa_example_smk_re
#' @description Calling `example("example_smk_re")` will run a random effects
#'   NMA model with the smoking cessation data, using the code in the Examples
#'   section below. The resulting `stan_nma` object `smk_fit_RE` will then be
#'   available in the global environment.
#'
#' @keywords examples
#' @template ex_smoking_network
#' @template ex_smoking_nma_re
#' @examples \dontshow{
#' if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) {
#'   assign("smk_net", smk_net, .GlobalEnv)
#'   assign("smk_fit_RE", smk_fit_RE, .GlobalEnv)
#' }
#' }
NULL

#' Example smoking UME NMA
#'
#' Smoking UME NMA for use in examples.
#'
#' @name example_smk_ume
#' @rdname aa_example_smk_ume
#' @description Calling `example("example_smk_ume")` will run an unrelated mean
#'   effects (inconsistency) NMA model with the smoking cessation data, using the
#'   code in the Examples section below. The resulting `stan_nma` object
#'   `smk_fit_RE_UME` will then be available in the global environment.
#'
#' @keywords examples
#' @template ex_smoking_network
#' @template ex_smoking_nma_re_ume
#' @examples \dontshow{
#' if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) {
#'   assign("smk_net", smk_net, .GlobalEnv)
#'   assign("smk_fit_RE_UME", smk_fit_RE_UME, .GlobalEnv)
#' }
#' }
NULL

#' Example plaque psoriasis ML-NMR
#'
#' Plaque psoriasis ML-NMR for use in examples.
#'
#' @name example_pso_mlnmr
#' @rdname aa_example_pso_mlnmr
#' @description Calling `example("example_pso_mlnmr")` will run a ML-NMR model
#'   with the plaque psoriasis IPD and AgD, using the code in the Examples
#'   section below. The resulting `stan_nma` object `pso_fit` will then be
#'   available in the global environment.
#'
#' @keywords examples
#' @template ex_plaque_psoriasis_network
#' @template ex_plaque_psoriasis_integration
#' @template ex_plaque_psoriasis_mlnmr
#' @examples \dontshow{
#' if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) {
#'   assign("pso_net", pso_net, .GlobalEnv)
#'   assign("pso_fit", pso_fit, .GlobalEnv)
#' }
#' }
NULL
