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
    male = as.numeric(male),
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                         trtn == 4 ~ "TNFa blocker",
                         trtn == 7 ~ "IL-12/23 blocker"),
    # Check complete cases for covariates of interest
    is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  ) %>%
  arrange(studyc, trtn)

test_that("drops IPD NA rows with warning by default", {
  ipd_net <- set_ipd(pso_ipd, study = studyc, trt = trtc, r = pasi75)

  expect_error(compare_populations(ipd_net), "provide `covariates` to compare on")
  expect_warning(comp <- compare_populations(ipd_net,
                                             covariates = c("durnpso", "prevsys", "psa", "bsa", "weight", "male", "age")),
                 "Removed 4 observations with missing covariate values from IPD")
  expect_s3_class(comp, "pop_comp")
  expect_equal(comp$method, "propensity")
})

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

pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male))

test_that("covariates must be listed in data", {
  m <- "Cannot compare requested covariates missing integration points"
  expect_error(compare_populations(pso_net, method = "propensity", covariates = "height"), m)
  expect_error(compare_populations(pso_net, method = "propensity", covariates = "a"), m)
  expect_error(compare_populations(pso_net, method = "propensity", covariates = c("age", "height")), m)
  expect_error(compare_populations(pso_net, method = "propensity", covariates = 1), m)
})

test_that("method argument checks", {
  expect_error(compare_populations(pso_net, method = 1), "must be a character vector")
  expect_error(compare_populations(pso_net, method = "a"), "must be one of")
})

test_that("output class", {
  out <- compare_populations(pso_net)

  expect_s3_class(out,"pop_comp")
  expect_equal(out$method, "propensity")
  expect_equal(unique(out$components$component), 1L)

  expect_equal(colnames(out$comparison_matrix), levels(pso_net$studies))
  expect_equal(rownames(out$comparison_matrix), levels(pso_net$studies))
  expect_equal(nrow(out$summary), choose(nlevels(pso_net$studies), 2))
})

test_that("print outputs", {
  comp <- compare_populations(pso_net)
  expect_output(print(comp),
                "Compared populations using propensity score overlap, based on the following covariates: durnpso, prevsys, bsa, weight, psa, age and male")
  expect_output(print(comp, n = 1),
                "UNCOVER-2 vs\\. UNCOVER-3 +2558 +2528\\.59 +98\\.85")
  expect_output(print(comp, n = 1, simplify = FALSE),
                "UNCOVER-2 vs\\. UNCOVER-3 +2558 +2528\\.59 +98\\.85")
})

# Disconnected network
pso_net <- combine_network(
  set_ipd(pso_ipd %>% filter(!trtc %in% c("SEC_150", "SEC_300", "UST")),
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass, allow_single_arm = TRUE),
  set_agd_arm(pso_agd %>% filter(trtc %in% c("SEC_150", "SEC_300", "UST")),
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass)
)

pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male))

test_that("output class", {
  out <- compare_populations(pso_net)

  expect_s3_class(out,"pop_comp")
  expect_equal(out$method, "propensity")
  expect_equal(unique(out$components$component), c(1L, 2L))

  expect_equal(colnames(out$comparison_matrix), levels(pso_net$studies))
  expect_equal(rownames(out$comparison_matrix), levels(pso_net$studies))
  expect_equal(nrow(out$summary), choose(nlevels(pso_net$studies), 2))
})

test_that("print outputs", {
  comp <- compare_populations(pso_net)
  expect_output(print(comp),
                "Compared populations using propensity score overlap, based on the following covariates: durnpso, prevsys, bsa, weight, psa, age and male")
  expect_output(print(comp, n = 1, simplify = TRUE),
                "UNCOVER-2 vs\\. CLEAR +1888 +1601\\.31 +84\\.81")
  expect_output(print(comp, n = 1, simplify = TRUE),
                "Subnetwork 2 vs\\. 1")
  expect_output(print(comp, n = 1, simplify = FALSE),
                "UNCOVER-2 vs\\. UNCOVER-3 +2558 +2528\\.59 +98\\.85")
})
