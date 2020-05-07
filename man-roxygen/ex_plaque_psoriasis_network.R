#' @examples
#' # Set up plaque psoriasis network combining IPD and AgD
#' library(dplyr)
#' pso_ipd <- filter(plaque_psoriasis_ipd,
#'                   studyc %in% c("UNCOVER-1", "UNCOVER-2", "UNCOVER-3"))
#'
#' pso_agd <- filter(plaque_psoriasis_agd,
#'                   studyc == "FIXTURE")
#'
#' head(pso_ipd)
#' head(pso_agd)
#'
#' pso_ipd <- pso_ipd %>%
#'   mutate(# Variable transformations
#'     bsa = bsa / 100,
#'     prevsys = as.numeric(prevsys),
#'     psa = as.numeric(psa),
#'     weight = weight / 10,
#'     durnpso = durnpso / 10,
#'     # Treatment classes
#'     trtclass = case_when(trtn == 1 ~ "Placebo",
#'                          trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
#'                          trtn == 4 ~ "TNFa blocker"),
#'     # Check complete cases for covariates of interest
#'     complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
#'   )
#'
#' pso_agd <- pso_agd %>%
#'   mutate(
#'     # Variable transformations
#'     bsa_mean = bsa_mean / 100,
#'     bsa_sd = bsa_sd / 100,
#'     prevsys = prevsys / 100,
#'     psa = psa / 100,
#'     weight_mean = weight_mean / 10,
#'     weight_sd = weight_sd / 10,
#'     durnpso_mean = durnpso_mean / 10,
#'     durnpso_sd = durnpso_sd / 10,
#'     # Treatment classes
#'     trtclass = case_when(trtn == 1 ~ "Placebo",
#'                          trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
#'                          trtn == 4 ~ "TNFa blocker")
#'   )
#'
#' # Exclude small number of individuals with missing covariates
#' pso_ipd <- filter(pso_ipd, complete)
#'
#' pso_net <- combine_network(
#'   set_ipd(pso_ipd,
#'           study = studyc,
#'           trt = trtc,
#'           r = pasi75,
#'           trt_class = trtclass),
#'   set_agd_arm(pso_agd,
#'               study = studyc,
#'               trt = trtc,
#'               r = pasi75_r,
#'               n = pasi75_n,
#'               trt_class = trtclass)
#' )
#'
#' # Print network details
#' pso_net
#'
