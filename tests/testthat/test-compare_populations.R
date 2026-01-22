library(multinma)
library(dplyr)

test_that("compare_populations() wants a nma_data object", {
  m <- "Expecting an `nma_data` object"
  expect_error(compare_populations(network = 1), m)
  expect_error(compare_populations(network = list()), m)
  expect_error(compare_populations(network = NULL), m)
})

pso_ipd <- plaque_psoriasis_ipd %>%
  mutate(
    # Variable transformations
    bsa = bsa / 100,
    weight = weight / 10,
    durnpso = durnpso / 10,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                         trtn == 4 ~ "TNFa blocker",
                         trtn == 7 ~ "IL-12/23 blocker"),
    # Check complete cases for covariates of interest
    is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  ) %>%
  arrange(studyc, trtn)

# AgD studies
pso_agd <- plaque_psoriasis_agd %>%
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    prevsys = prevsys / 100,
    psa = psa / 100,
    male = male / 100,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                         trtn == 4 ~ "TNFa blocker",
                         trtn == 7 ~ "IL-12/23 blocker")
  ) %>%
  arrange(studyc, trtn)


# Missing Data
pso_ipd %>%
  group_by(studyc) %>%
  summarise(n_total = n(),
            n_missing = sum(!is_complete),
            pct_missing = mean(!is_complete) * 100)

pso_ipd <- filter(pso_ipd, is_complete)

# Creating the FULL network
pso_net <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass)
)

test_that("compare_populations() covariates argument must be one listed", {
  m <- "The network object does not contain integration code"
  expect_error(compare_populations(pso_net), m)
})

pso_net_int <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male),
                           n_int = 64)

test_that("compare_populations() covariates must be listed in data", {
  m <- "The following covariates are NOT defined in the network's integration code"
  expect_error(compare_populations(pso_net_int, covariates = "height"), m)
  expect_error(compare_populations(pso_net_int, covariates = "a"), m)
  expect_error(compare_populations(pso_net_int, covariates = c("age", "height")), m)
  expect_error(compare_populations(pso_net_int, covariates = 1), m)
})

test_that("a method is stated in compare_populations()", {
  m <- "Argument `method` must be either 'euclidean' or 'propensity'."
  expect_error(compare_populations(pso_net_int, method = 1), m)
  expect_error(compare_populations(pso_net_int, method = "a"), m)
  expect_error(compare_populations(pso_net_int, method = euclidean), m)
  expect_error(compare_populations(pso_net_int, method = c("euclidean", "propensity")), m)
})

#-------------------------------------------------
#---Cross validation tests
#-------------------------------------------------

fit <- nma(pso_net_int,
           trt_effects = "fixed",
           link = "probit",
           regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10),
           prior_aux = flat(),
           QR = TRUE,
           init_r = 0.5,
           chains = 2, iter = 20, warmup = 10)


test_that("cross_validation() rejects incorrect input types", {
  m <- "Input must be a 'stan_nma' object"
  expect_error(cross_validation(NULL), m)
  expect_error(cross_validation(list(a = 1)), m)
  expect_error(cross_validation(pso_net), m)
})


