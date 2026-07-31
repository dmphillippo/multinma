library(multinma)
library(dplyr)

test_that("compare_populations() wants a nma_data object", {
  m <- "Expecting an `nma_data` object"
  expect_error(compare_populations(network = 1, method = "propensity"), m)
  expect_error(compare_populations(network = list(), method = "propensity"), m)
  expect_error(compare_populations(network = NULL, method = "propensity"), m)
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

test_that("type = 'propensity' needs integration points for AgD", {
  m <- "Integration points must be present"
  expect_error(compare_populations(pso_net, method = "propensity"), m)
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
  m <- "Cannot compare requested covariates missing integration points"
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = "height"), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = "a"), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = c("age", "height")), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = 1), m)
})

test_that("a method is stated in compare_populations()", {
  expect_error(compare_populations(pso_net_int, method = 1), "must be a character vector")
  expect_error(compare_populations(pso_net_int, method = "a"), "must be one of")
})
