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

  expect_error(compare_populations(ipd_net), "Provide `covariates` to compare on")
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
    bsa = bsa_mean, # add in means again for euclidean method
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    weight = weight_mean,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    durnpso = durnpso_mean,
    age = age_mean,
    prevsys = prevsys / 100,
    prevsys_sd = sqrt(prevsys * (1 - prevsys)),
    psa = psa / 100,
    psa_sd = sqrt(psa * (1 - psa)),
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

test_that("type = 'euclidean'", {
  expect_error(compare_populations(pso_net, method = "euclidean"),
               'Provide `covariates` to compare on when method = "euclidean"')

  # _sd column missing
  expect_error(compare_populations(pso_net, method = "euclidean",
                                   covariates = c("weight", "male")),
               "Standard deviation columns not found in AgD: male_sd")

  # with sd column
  comp <- compare_populations(pso_net, method = "euclidean",
            covariates = c("durnpso", "prevsys", "psa", "bsa", "weight", "age"))

  expect_s3_class(comp, "pop_comp")
  expect_equal(comp$method, "euclidean")
  expect_output(print(comp, n = 1),
                "CLEAR vs\\. ERASURE +1403 +0\\.17")
  expect_output(print(comp, n = 1, order = "increasing"),
                "IXORA-S vs\\. JUNCTURE +441 +0\\.92")
})

pso_net_int <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male))

test_that("covariates must be listed in data", {
  m <- "Cannot compare requested covariates missing integration points"
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = "height"), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = "a"), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = c("age", "height")), m)
  expect_error(compare_populations(pso_net_int, method = "propensity", covariates = 1), m)

  expect_error(compare_populations(pso_net_int, method = "euclidean", covariates = "height"), m)
  expect_error(compare_populations(pso_net_int, method = "euclidean", covariates = "a"), m)
  expect_error(compare_populations(pso_net_int, method = "euclidean", covariates = c("age", "height")), m)
  expect_error(compare_populations(pso_net_int, method = "euclidean", covariates = 1), m)
})

test_that("method argument checks", {
  expect_error(compare_populations(pso_net_int, method = 1), "must be a character vector")
  expect_error(compare_populations(pso_net_int, method = "a"), "must be one of")
})

test_that("output class", {
  out <- compare_populations(pso_net_int)

  expect_s3_class(out,"pop_comp")
  expect_equal(out$method, "propensity")
  expect_equal(unique(out$components$component), 1L)

  expect_equal(colnames(out$comparison_matrix), levels(pso_net_int$studies))
  expect_equal(rownames(out$comparison_matrix), levels(pso_net_int$studies))
  expect_equal(nrow(out$summary), choose(nlevels(pso_net_int$studies), 2))

  out <- compare_populations(pso_net_int, method = "euclidean")

  expect_s3_class(out,"pop_comp")
  expect_equal(out$method, "euclidean")
  expect_equal(unique(out$components$component), 1L)

  expect_equal(colnames(out$comparison_matrix), levels(pso_net_int$studies))
  expect_equal(rownames(out$comparison_matrix), levels(pso_net_int$studies))
  expect_equal(nrow(out$summary), choose(nlevels(pso_net_int$studies), 2))
})

test_that("print outputs", {
  comp <- compare_populations(pso_net_int)
  expect_output(print(comp),
                "Compared populations using propensity score overlap, based on the following covariates: durnpso, prevsys, bsa, weight, psa, age and male")
  expect_output(print(comp, n = 1),
                "UNCOVER-2 vs\\. UNCOVER-3 +2558 +2534\\.98 +99\\.1")
  expect_output(print(comp, n = 1, simplify = FALSE),
                "UNCOVER-2 vs\\. UNCOVER-3 +2558 +2534\\.98 +99\\.1")
})

test_that("type = 'euclidean' with integration points", {
  comp <- compare_populations(pso_net_int, method = "euclidean")

  expect_s3_class(comp, "pop_comp")
  expect_equal(comp$method, "euclidean")
  expect_output(print(comp, n = 1),
                "CLEAR vs\\. ERASURE +1403 +0\\.19")
  expect_output(print(comp, n = 1, order = "increasing"),
                "IXORA-S vs\\. JUNCTURE +441 +0\\.93")

  # equals non-integration version
  comp_int <- compare_populations(pso_net_int, method = "euclidean",
                              covariates = c("durnpso", "prevsys", "psa", "bsa", "weight", "age"))

  comp_stat <- compare_populations(pso_net, method = "euclidean",
                              covariates = c("durnpso", "prevsys", "psa", "bsa", "weight", "age"))

  expect_equal(comp_int$comparison_matrix, comp_stat$comparison_matrix, tol = 0.05)
})

# Disconnected network
pso_net_disc <- combine_network(
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

pso_net_disc <- add_integration(pso_net_disc,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male))

test_that("output class", {
  out <- compare_populations(pso_net_disc)

  expect_s3_class(out,"pop_comp")
  expect_equal(out$method, "propensity")
  expect_equal(unique(out$components$component), c(1L, 2L))

  expect_equal(colnames(out$comparison_matrix), levels(pso_net_disc$studies))
  expect_equal(rownames(out$comparison_matrix), levels(pso_net_disc$studies))
  expect_equal(nrow(out$summary), choose(nlevels(pso_net_disc$studies), 2))
})

test_that("print outputs", {
  comp <- compare_populations(pso_net_disc)
  expect_output(print(comp),
                "Compared populations using propensity score overlap, based on the following covariates: durnpso, prevsys, bsa, weight, psa, age and male")
  expect_output(print(comp, n = 1, simplify = TRUE),
                "ERASURE vs\\. UNCOVER-3 +1827 +1781\\.52 +97\\.51")
  expect_output(print(comp, n = 1, simplify = TRUE),
                "Subnetwork 2 vs\\. 1")
  expect_output(print(comp, n = 1, simplify = FALSE),
                "CLEAR vs\\. ERASURE +1157 +1151\\.69 +99\\.54")
})

