
smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Only test gradients, no sampling
smk_fit_RE <- nma(smk_net,
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = normal(scale = 5),
                  test_grad = TRUE)

test_that("all_contrasts argument", {
  m <- "should be TRUE or FALSE"
  expect_error(marginal_effects(smk_fit_RE, all_contrasts = "a"), m)
  expect_error(marginal_effects(smk_fit_RE, all_contrasts = 1), m)
  expect_error(marginal_effects(smk_fit_RE, all_contrasts = list()), m)
  expect_error(marginal_effects(smk_fit_RE, all_contrasts = NA), m)
  expect_error(marginal_effects(smk_fit_RE, all_contrasts = NULL), m)
})

test_that("trt_ref argument", {
  m <- "does not match a treatment in the network"
  expect_error(marginal_effects(smk_fit_RE, trt_ref = "a"), m)
  expect_error(marginal_effects(smk_fit_RE, trt_ref = 1), m)
  expect_error(marginal_effects(smk_fit_RE, trt_ref = list("a")), m)
  expect_error(marginal_effects(smk_fit_RE, trt_ref = NA), m)

  # expect_warning(marginal_effects(smk_fit_RE, all_contrasts = TRUE, trt_ref = "Self-help"), "Ignoring `trt_ref`")
})

test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(marginal_effects(smk_fit_RE, probs = "a"), m)
  expect_error(marginal_effects(smk_fit_RE, probs = -1), m)
  expect_error(marginal_effects(smk_fit_RE, probs = 1.5), m)
  expect_error(marginal_effects(smk_fit_RE, probs = Inf), m)
  expect_error(marginal_effects(smk_fit_RE, probs = list()), m)
  expect_error(marginal_effects(smk_fit_RE, probs = NA), m)
  expect_error(marginal_effects(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(marginal_effects(smk_fit_RE, summary = "a"), m)
  expect_error(marginal_effects(smk_fit_RE, summary = 1), m)
  expect_error(marginal_effects(smk_fit_RE, summary = list()), m)
  expect_error(marginal_effects(smk_fit_RE, summary = NA), m)
  expect_error(marginal_effects(smk_fit_RE, summary = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(marginal_effects(smk_fit_RE, newdata = "a"), m)
  expect_error(marginal_effects(smk_fit_RE, newdata = 1), m)
  expect_error(marginal_effects(smk_fit_RE, newdata = list()), m)
  expect_error(marginal_effects(smk_fit_RE, newdata = NA), m)
})


skip_on_cran()  # Reduce CRAN check time


# Binary AgD only --------------------------------------------------------------

# Only small number of samples to test
smk_fit_RE <- suppressWarnings(nma(smk_net,
                                   trt_effects = "random",
                                   prior_intercept = normal(scale = 100),
                                   prior_trt = normal(scale = 100),
                                   prior_het = normal(scale = 5),
                                   iter = 10))

test_that("Binary AgD: .trta, .trtb columns are correct", {
  re1.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re1.1$.trtb, "]"),
                   re1.1$parameter)
  expect_identical(as.character(re1.1$.trta), rep("No intervention", times = nrow(re1.1)))

  re1.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, mtype = "ratio"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re1.2$.trtb, "]"),
                   re1.2$parameter)
  expect_identical(as.character(re1.2$.trta), rep("No intervention", times = nrow(re1.2)))

  re1.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, mtype = "link"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re1.3$.trtb, "]"),
                   re1.3$parameter)
  expect_identical(as.character(re1.3$.trta), rep("No intervention", times = nrow(re1.3)))


  re2.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, all_contrasts = TRUE, mtype = "difference"))
  expect_identical(paste0("marg[", re2.1$.study, ": ", re2.1$.trtb, " vs. ", re2.1$.trta, "]"),
                   re2.1$parameter)

  re2.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, all_contrasts = TRUE, mtype = "ratio"))
  expect_identical(paste0("marg[", re2.1$.study, ": ", re2.2$.trtb, " vs. ", re2.2$.trta, "]"),
                   re2.2$parameter)

  re2.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, all_contrasts = TRUE, mtype = "link"))
  expect_identical(paste0("marg[", re2.1$.study, ": ", re2.3$.trtb, " vs. ", re2.3$.trta, "]"),
                   re2.3$parameter)


  re3.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, trt_ref = "Self-help", mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re3.1$.trtb, "]"),
                   re3.1$parameter)
  expect_identical(as.character(re3.1$.trta), rep("Self-help", times = nrow(re3.1)))

  re3.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, trt_ref = "Self-help", mtype = "ratio"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re3.2$.trtb, "]"),
                   re3.2$parameter)
  expect_identical(as.character(re3.2$.trta), rep("Self-help", times = nrow(re3.2)))

  re3.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, trt_ref = "Self-help", mtype = "link"))
  expect_identical(paste0("marg[", re3.1$.study, ": ", re3.3$.trtb, "]"),
                   re3.3$parameter)
  expect_identical(as.character(re3.3$.trta), rep("Self-help", times = nrow(re3.3)))
})


test_that("Binary AgD: .trta, .trtb columns are correct in new populations", {
  re1.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.trtb, "]"),
                   re1.1$parameter)
  expect_identical(as.character(re1.1$.trta), rep("No intervention", times = nrow(re1.1)))

  re1.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), mtype = "ratio"))
  expect_identical(paste0("marg[", re1.2$.trtb, "]"),
                   re1.2$parameter)
  expect_identical(as.character(re1.2$.trta), rep("No intervention", times = nrow(re1.2)))

  re1.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), mtype = "link"))
  expect_identical(paste0("marg[", re1.3$.trtb, "]"),
                   re1.3$parameter)
  expect_identical(as.character(re1.3$.trta), rep("No intervention", times = nrow(re1.3)))


  re2.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "difference"))
  expect_identical(paste0("marg[", re2.1$.trtb, " vs. ", re2.1$.trta, "]"),
                   re2.1$parameter)

  re2.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "ratio"))
  expect_identical(paste0("marg[", re2.2$.trtb, " vs. ", re2.2$.trta, "]"),
                   re2.2$parameter)

  re2.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "link"))
  expect_identical(paste0("marg[", re2.3$.trtb, " vs. ", re2.3$.trta, "]"),
                   re2.3$parameter)


  re3.1 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), trt_ref = "Self-help", mtype = "difference"))
  expect_identical(paste0("marg[", re3.1$.trtb, "]"),
                   re3.1$parameter)
  expect_identical(as.character(re3.1$.trta), rep("Self-help", times = nrow(re3.1)))

  re3.2 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), trt_ref = "Self-help", mtype = "ratio"))
  expect_identical(paste0("marg[", re3.2$.trtb, "]"),
                   re3.2$parameter)
  expect_identical(as.character(re3.2$.trta), rep("Self-help", times = nrow(re3.2)))

  re3.3 <- tibble::as_tibble(marginal_effects(smk_fit_RE, baseline = distr(qnorm, 0, 0.1), trt_ref = "Self-help", mtype = "link"))
  expect_identical(paste0("marg[", re3.3$.trtb, "]"),
                   re3.3$parameter)
  expect_identical(as.character(re3.3$.trta), rep("Self-help", times = nrow(re3.3)))
})

blocker_net <- set_agd_arm(blocker, studyn, trtc, r = r, n = n)
blocker_fit <- suppressWarnings(nma(blocker_net,
                                    prior_intercept = normal(scale = 100),
                                    prior_trt = normal(scale = 100),
                                    iter = 10))

test_that("all_contrasts works with only two treatments", {
  expect_equal(
    dplyr::select(as.data.frame(marginal_effects(blocker_fit)), -"parameter"),
    dplyr::select(as.data.frame(marginal_effects(blocker_fit, all_contrasts = TRUE)) ,-"parameter"),
    check.attributes = FALSE
  )
  expect_identical(
    unname(as.array(marginal_effects(blocker_fit))),
    unname(as.array(marginal_effects(blocker_fit, all_contrasts = TRUE)))
  )
})


# Binary with IPD regression ---------------------------------------------------

library(dplyr)

pso_ipd <- mutate(plaque_psoriasis_ipd,
    bsa = bsa / 100,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    weight = weight / 10,
    durnpso = durnpso / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker"),
    # Check complete cases for covariates of interest
    complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  )

pso_agd <- mutate(plaque_psoriasis_agd,
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    prevsys = prevsys / 100,
    psa = psa / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker")
  )

pso_net <- combine_network(
  set_ipd(filter(pso_ipd, complete),
          studyc, trtc,
          r = pasi75),
  set_agd_arm(pso_agd,
              study = studyc, trt = trtc,
              r = pasi75_r, n = pasi75_n))

pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           n_int = 4)

# Only small number of samples to test
pso_fit <- suppressWarnings(nma(pso_net,
                                trt_effects = "fixed",
                                regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
                                prior_intercept = normal(scale = 10),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                init_r = 0.1,
                                iter = 10))

pso_new <- tibble(
  bsa_mean = 0.6,
  bsa_sd = 0.3,
  prevsys = 0.1,
  psa = 0.2,
  weight_mean = 10,
  weight_sd = 1,
  durnpso_mean = 3,
  durnpso_sd = 1,
  studyc = "NEW"
)

pso_new <- add_integration(pso_new,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           cor = pso_net$int_cor,
                           n_int = 4)


test_that("newdata validation", {
  expect_error(marginal_effects(pso_fit, newdata = dplyr::select(pso_new, -".int_durnpso"), baseline = distr(qnorm, 0, 0.1)),
               'Regression variable "durnpso" not found in `newdata`')
  expect_error(marginal_effects(pso_fit, newdata = dplyr::select(pso_new, -".int_durnpso", -".int_weight"), baseline = distr(qnorm, 0, 0.1)),
               'Regression variables "durnpso" and "weight" not found in `newdata`')
})

test_that("baseline required when newdata specified", {
  m <- "Specify both `newdata` and `baseline`, or neither"
  expect_error(marginal_effects(pso_fit, newdata = pso_new), m)
  expect_error(marginal_effects(pso_fit, baseline = distr(qnorm, 0, 0.1)), m)
})

test_that("Binary IPD: .study, .trta, .trtb columns are correct", {
  re1 <- tibble::as_tibble(marginal_effects(pso_fit))
  expect_identical(paste0("marg[", re1$.study, ": ", re1$.trtb, "]"),
                   re1$parameter)
  expect_identical(as.character(re1$.trta), rep(levels(pso_net$treatments)[1], times = nrow(re1)))

  re2 <- tibble::as_tibble(marginal_effects(pso_fit, all_contrasts = TRUE))
  expect_identical(paste0("marg[", re2$.study, ": ", re2$.trtb, " vs. ", re2$.trta, "]"),
                   re2$parameter)

  re3 <- tibble::as_tibble(marginal_effects(pso_fit, trt_ref = "ETN"))
  expect_identical(paste0("marg[", re3$.study, ": ", re3$.trtb, "]"),
                   re3$parameter)
  expect_identical(as.character(re3$.trta), rep("ETN", times = nrow(re3)))

  re4 <- tibble::as_tibble(marginal_effects(pso_fit, baseline = distr(qnorm, 0, 0.1), newdata = pso_new, study = studyc))
  expect_identical(paste0("marg[", re4$.study, ": ", re4$.trtb, "]"),
                   re4$parameter)
  expect_identical(as.character(re4$.trta), rep(levels(pso_net$treatments)[1], times = nrow(re4)))

  re5 <- tibble::as_tibble(marginal_effects(pso_fit, all_contrasts = TRUE, baseline = distr(qnorm, 0, 0.1), newdata = pso_new, study = studyc))
  expect_identical(paste0("marg[", re5$.study, ": ", re5$.trtb, " vs. ", re5$.trta, "]"),
                   re5$parameter)

  re6 <- tibble::as_tibble(marginal_effects(pso_fit, trt_ref = "ETN", baseline = distr(qnorm, 0, 0.1), newdata = pso_new, study = studyc))
  expect_identical(paste0("marg[", re6$.study, ": ", re6$.trtb, "]"),
                   re6$parameter)
  expect_identical(as.character(re6$.trta), rep("ETN", times = nrow(re6)))
})


# Ordered categorical -----------------------------------------------------

ord_net <- set_agd_arm(plaque_psoriasis_agd,
                       studyc, trtc,
                       r = multi(r0 = pasi75_n,
                                 PASI75 = pasi75_r,
                                 PASI90 = pasi90_r,
                                 PASI100 = pasi100_r,
                                 type = "ordered", inclusive = TRUE))


# Only small number of samples to test
ord_fit <- suppressWarnings(nma(ord_net,
                                trt_effects = "random",
                                prior_intercept = normal(scale = 100),
                                prior_trt = normal(scale = 100),
                                prior_het = normal(scale = 5),
                                iter = 10))

test_that("Ordered categorical: .study, .trta, .trtb, .category columns are correct", {
  re1.1 <- tibble::as_tibble(marginal_effects(ord_fit, mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re1.1$.trtb, ", ", re1.1$.category, "]"),
                   re1.1$parameter)
  expect_identical(as.character(re1.1$.trta), rep("SEC_300", times = nrow(re1.1)))

  re1.2 <- tibble::as_tibble(marginal_effects(ord_fit, mtype = "ratio"))
  expect_identical(paste0("marg[", re1.2$.study, ": ", re1.2$.trtb, ", ", re1.2$.category, "]"),
                   re1.2$parameter)
  expect_identical(as.character(re1.2$.trta), rep("SEC_300", times = nrow(re1.2)))

  re1.3 <- tibble::as_tibble(marginal_effects(ord_fit, mtype = "link"))
  expect_identical(paste0("marg[", re1.3$.study, ": ", re1.3$.trtb, ", ", re1.3$.category, "]"),
                   re1.3$parameter)
  expect_identical(as.character(re1.3$.trta), rep("SEC_300", times = nrow(re1.3)))


  re2.1 <- tibble::as_tibble(marginal_effects(ord_fit, all_contrasts = TRUE, mtype = "difference"))
  expect_identical(paste0("marg[", re2.1$.study, ": ", re2.1$.trtb, " vs. ", re2.1$.trta, ", ", re2.1$.category, "]"),
                   re2.1$parameter)

  re2.2 <- tibble::as_tibble(marginal_effects(ord_fit, all_contrasts = TRUE, mtype = "ratio"))
  expect_identical(paste0("marg[", re2.2$.study, ": ", re2.2$.trtb, " vs. ", re2.2$.trta, ", ", re2.2$.category, "]"),
                   re2.2$parameter)

  re2.3 <- tibble::as_tibble(marginal_effects(ord_fit, all_contrasts = TRUE, mtype = "link"))
  expect_identical(paste0("marg[", re2.3$.study, ": ", re2.3$.trtb, " vs. ", re2.3$.trta, ", ", re2.3$.category, "]"),
                   re2.3$parameter)


  re3.1 <- tibble::as_tibble(marginal_effects(ord_fit, trt_ref = "PBO", mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re3.1$.trtb, ", ", re3.1$.category, "]"),
                   re3.1$parameter)
  expect_identical(as.character(re3.1$.trta), rep("PBO", times = nrow(re3.1)))

  re3.2 <- tibble::as_tibble(marginal_effects(ord_fit, trt_ref = "PBO", mtype = "ratio"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re3.2$.trtb, ", ", re3.2$.category, "]"),
                   re3.2$parameter)
  expect_identical(as.character(re3.2$.trta), rep("PBO", times = nrow(re3.2)))

  re3.3 <- tibble::as_tibble(marginal_effects(ord_fit, trt_ref = "PBO", mtype = "link"))
  expect_identical(paste0("marg[", re3.1$.study, ": ", re3.3$.trtb, ", ", re3.3$.category, "]"),
                   re3.3$parameter)
  expect_identical(as.character(re3.3$.trta), rep("PBO", times = nrow(re3.3)))
})


test_that("Ordered categorical: .trta, .trtb, .category columns are correct in new populations", {
  re1.1 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.trtb, ", ", re1.1$.category, "]"),
                   re1.1$parameter)
  expect_identical(as.character(re1.1$.trta), rep("SEC_300", times = nrow(re1.1)))

  re1.2 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), mtype = "ratio"))
  expect_identical(paste0("marg[", re1.2$.trtb, ", ", re1.2$.category, "]"),
                   re1.2$parameter)
  expect_identical(as.character(re1.2$.trta), rep("SEC_300", times = nrow(re1.2)))

  re1.3 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), mtype = "link"))
  expect_identical(paste0("marg[", re1.3$.trtb, ", ", re1.3$.category, "]"),
                   re1.3$parameter)
  expect_identical(as.character(re1.3$.trta), rep("SEC_300", times = nrow(re1.3)))


  re2.1 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "difference"))
  expect_identical(paste0("marg[", re2.1$.trtb, " vs. ", re2.1$.trta, ", ", re2.1$.category, "]"),
                   re2.1$parameter)

  re2.2 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "ratio"))
  expect_identical(paste0("marg[", re2.2$.trtb, " vs. ", re2.2$.trta, ", ", re2.2$.category, "]"),
                   re2.2$parameter)

  re2.3 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), all_contrasts = TRUE, mtype = "link"))
  expect_identical(paste0("marg[", re2.3$.trtb, " vs. ", re2.3$.trta, ", ", re2.3$.category, "]"),
                   re2.3$parameter)


  re3.1 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), trt_ref = "PBO", mtype = "difference"))
  expect_identical(paste0("marg[", re3.1$.trtb, ", ", re3.1$.category, "]"),
                   re3.1$parameter)
  expect_identical(as.character(re3.1$.trta), rep("PBO", times = nrow(re3.1)))

  re3.2 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), trt_ref = "PBO", mtype = "ratio"))
  expect_identical(paste0("marg[", re3.2$.trtb, ", ", re3.2$.category, "]"),
                   re3.2$parameter)
  expect_identical(as.character(re3.2$.trta), rep("PBO", times = nrow(re3.2)))

  re3.3 <- tibble::as_tibble(marginal_effects(ord_fit, baseline = distr(qnorm, 0, 0.1), trt_ref = "PBO", mtype = "link"))
  expect_identical(paste0("marg[", re3.3$.trtb, ", ", re3.3$.category, "]"),
                   re3.3$parameter)
  expect_identical(as.character(re3.3$.trta), rep("PBO", times = nrow(re3.3)))
})


# Survival ----------------------------------------------------------------

surv_net <- set_ipd(ndmm_agd, study, trt, Surv = Surv(eventtime, status))

# Only small number of samples to test
surv_fit <- suppressWarnings(nma(surv_net,
                                 likelihood = "weibull",
                                 prior_intercept = normal(scale = 100),
                                 prior_trt = normal(scale = 100),
                                 prior_aux = normal(scale = 5),
                                 iter = 10))

test_that("Survival: .study, .trta, .trtb columns are correct", {
  re1.1 <- tibble::as_tibble(marginal_effects(surv_fit, type = "median", mtype = "difference"))
  expect_identical(paste0("marg[", re1.1$.study, ": ", re1.1$.trtb, "]"),
                   re1.1$parameter)
  expect_identical(as.character(re1.1$.trta), rep("Pbo", times = nrow(re1.1)))

  re1.2 <- tibble::as_tibble(marginal_effects(surv_fit, type = "hazard", mtype = "ratio")) %>%
    group_by(.study, .trtb, .trta) %>%
    mutate(.id = 1:n())
  expect_identical(paste0("marg[", re1.2$.study, ": ", re1.2$.trtb, ", ", re1.2$.id, "]"),
                   re1.2$parameter)
  expect_identical(as.character(re1.2$.trta), rep("Pbo", times = nrow(re1.2)))

  re1.3 <- tibble::as_tibble(marginal_effects(surv_fit, type = "quantile", mtype = "difference"))
  expect_identical(paste0("marg[", re1.3$.study, ": ", re1.3$.trtb, ", ", re1.3$.quantile, "]"),
                   re1.3$parameter)
  expect_identical(as.character(re1.3$.trta), rep("Pbo", times = nrow(re1.3)))


  re2.1 <- tibble::as_tibble(marginal_effects(surv_fit, all_contrasts = TRUE, type = "median", mtype = "difference"))
  expect_identical(paste0("marg[", re2.1$.study, ": ", re2.1$.trtb, " vs. ", re2.1$.trta, "]"),
                   re2.1$parameter)

  re2.2 <- tibble::as_tibble(marginal_effects(surv_fit, all_contrasts = TRUE, type = "hazard", mtype = "ratio")) %>%
    group_by(.study, .trtb, .trta) %>%
    mutate(.id = 1:n())
  expect_identical(paste0("marg[", re2.2$.study, ": ", re2.2$.trtb, " vs. ", re2.2$.trta, ", ", re2.2$.id, "]"),
                   re2.2$parameter)

  re2.3 <- tibble::as_tibble(marginal_effects(surv_fit, all_contrasts = TRUE, type = "quantile", mtype = "difference"))
  expect_identical(paste0("marg[", re2.3$.study, ": ", re2.3$.trtb, " vs. ", re2.3$.trta, ", ", re2.3$.quantile, "]"),
                   re2.3$parameter)


  re3.1 <- tibble::as_tibble(marginal_effects(surv_fit, trt_ref = "Len", type = "median", mtype = "difference"))
  expect_identical(paste0("marg[", re3.1$.study, ": ", re3.1$.trtb, "]"),
                   re3.1$parameter)
  expect_identical(as.character(re3.1$.trta), rep("Len", times = nrow(re3.1)))

  re3.2 <- tibble::as_tibble(marginal_effects(surv_fit, trt_ref = "Len", type = "hazard", mtype = "ratio")) %>%
    group_by(.study, .trtb, .trta) %>%
    mutate(.id = 1:n())
  expect_identical(paste0("marg[", re3.2$.study, ": ", re3.2$.trtb, ", ", re3.2$.id, "]"),
                   re3.2$parameter)
  expect_identical(as.character(re3.2$.trta), rep("Len", times = nrow(re3.2)))

  re3.3 <- tibble::as_tibble(marginal_effects(surv_fit, trt_ref = "Len", type = "quantile", mtype = "difference"))
  expect_identical(paste0("marg[", re3.3$.study, ": ", re3.3$.trtb, ", ", re3.3$.quantile, "]"),
                   re3.3$parameter)
  expect_identical(as.character(re3.3$.trta), rep("Len", times = nrow(re3.3)))
})
