
smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Only test gradients, no sampling
smk_fit_RE <- suppressWarnings(nma(smk_net,
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = normal(scale = 5),
                  iter = 10))

test_that("baseline argument", {
  m <- "should be specified using distr"
  expect_error(predict(smk_fit_RE, baseline = 1), m)
  expect_error(predict(smk_fit_RE, baseline = list("a")), m)
  expect_error(predict(smk_fit_RE, baseline = NA), m)

  expect_error(predict(smk_fit_RE, baseline = "a"), "`baseline` must match the name of an IPD or AgD \\(arm-based\\) study in the network")
})

test_that("trt_ref argument", {
  m <- "does not match a treatment in the network"
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = "a"), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = 1), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = list("a")), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = NA), m)
})

test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(predict(smk_fit_RE, probs = "a"), m)
  expect_error(predict(smk_fit_RE, probs = -1), m)
  expect_error(predict(smk_fit_RE, probs = 1.5), m)
  expect_error(predict(smk_fit_RE, probs = Inf), m)
  expect_error(predict(smk_fit_RE, probs = list()), m)
  expect_error(predict(smk_fit_RE, probs = NA), m)
  expect_error(predict(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(predict(smk_fit_RE, summary = "a"), m)
  expect_error(predict(smk_fit_RE, summary = 1), m)
  expect_error(predict(smk_fit_RE, summary = list()), m)
  expect_error(predict(smk_fit_RE, summary = NA), m)
  expect_error(predict(smk_fit_RE, summary = NULL), m)
})

test_that("expand argument", {
  m <- "should be TRUE or FALSE"
  expect_error(predict(smk_fit_RE, expand = "a"), m)
  expect_error(predict(smk_fit_RE, expand = 1), m)
  expect_error(predict(smk_fit_RE, expand = list()), m)
  expect_error(predict(smk_fit_RE, expand = NA), m)
  expect_error(predict(smk_fit_RE, expand = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(predict(smk_fit_RE, newdata = "a"), m)
  expect_error(predict(smk_fit_RE, newdata = 1), m)
  expect_error(predict(smk_fit_RE, newdata = list()), m)
  expect_error(predict(smk_fit_RE, newdata = NA), m)
})

test_that("type argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, type = "a"), m)
  expect_error(predict(smk_fit_RE, type = "lin"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, type = 1), m2)
  expect_error(predict(smk_fit_RE, type = list("a")), m2)
  expect_error(predict(smk_fit_RE, type = NA), m2)
})

test_that("level argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, level = "a"), m)
  expect_error(predict(smk_fit_RE, level = "agg"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, level = 1), m2)
  expect_error(predict(smk_fit_RE, level = list("a")), m2)
  expect_error(predict(smk_fit_RE, level = NA), m2)

  expect_error(predict(smk_fit_RE, level = "individual"),
               "Cannot produce individual predictions without a regression model.")
})

test_that("baseline_type argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, baseline_type = "a"), m)
  expect_error(predict(smk_fit_RE, baseline_type = "lin"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, baseline_type = 1), m2)
  expect_error(predict(smk_fit_RE, baseline_type = list("a")), m2)
  expect_error(predict(smk_fit_RE, baseline_type = NA), m2)
})

test_that("baseline_level argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, baseline_level = "a"), m)
  expect_error(predict(smk_fit_RE, baseline_level = "agg"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, baseline_level = 1), m2)
  expect_error(predict(smk_fit_RE, baseline_level = list("a")), m2)
  expect_error(predict(smk_fit_RE, baseline_level = NA), m2)
})


skip_on_cran()  # Reduce CRAN check time

test_that(".study, .trt columns are correct", {
  pred1 <- tibble::as_tibble(predict(smk_fit_RE))
  expect_identical(paste0("pred[", pred1$.study, ": ", pred1$.trt, "]"),
                   pred1$parameter)

  pred2 <- tibble::as_tibble(predict(smk_fit_RE, baseline = distr(qnorm, 0, 1)))
  expect_identical(paste0("pred[", pred2$.trt, "]"),
                   pred2$parameter)
})

test_that("expand TRUE and FALSE", {
  pred1 <- tibble::as_tibble(predict(smk_fit_RE, expand = TRUE))
  expect_identical(pred1$parameter,
                   paste0("pred[",
                          tidyr::expand_grid(s = smk_net$studies, t = smk_net$treatments) %>%
                            dplyr::mutate(l = paste(s, t, sep = ": ")) %>%
                            dplyr::pull(l), "]"))

  pred2 <- tibble::as_tibble(predict(smk_fit_RE, expand = FALSE))
  expect_identical(pred2$parameter,
                   paste0("pred[",
                          dplyr::arrange(smk_net$agd_arm, .study, .trt) %>%
                            dplyr::mutate(l = paste(.study, .trt, sep = ": ")) %>%
                            dplyr::pull(l), "]"))

  pred3 <- tibble::as_tibble(predict(smk_fit_RE, expand = FALSE, baseline = "1"))
  expect_identical(pred3$parameter,
                   paste0("pred[",
                          dplyr::filter(smk_net$agd_arm, .study == "1")$.trt, "]"))

  expect_error(predict(smk_fit_RE, expand = FALSE, baseline = distr(qnorm, 0, 1)),
               "Provide `newdata` with a `\\.trt` column")

  pred4 <- tibble::as_tibble(predict(smk_fit_RE, expand = FALSE,
                                     baseline = distr(qnorm, 0, 1),
                                     newdata = data.frame(.trt = c("No intervention", "Self-help"))))
  expect_identical(pred4$parameter,
                   paste0("pred[", c("No intervention", "Self-help"), "]"))

  expect_error(predict(smk_fit_RE, expand = FALSE,
                       baseline = distr(qnorm, 0, 1),
                       newdata = data.frame(.trt = c("a", "b"))),
               "Treatments in `newdata` do not match those in the network")
})

pso_net <- set_ipd(plaque_psoriasis_ipd[complete.cases(plaque_psoriasis_ipd), ],
                   studyc, trtc,
                   r = pasi75)

# Only small number of samples to test
pso_fit <- suppressWarnings(nma(pso_net,
               trt_effects = "fixed",
               regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10),
               init_r = 0.1,
               iter = 10))

pso_new <- data.frame(durnpso = 10, prevsys = TRUE, bsa = 20, weight = 75, psa = FALSE, study = c("One", "Two"))

test_that("baseline and newdata for regression models", {
  m <- "Specify both `newdata` and `baseline`, or neither"
  expect_error(predict(pso_fit, newdata = pso_new), m)
  expect_error(predict(pso_fit, baseline = distr(qnorm, 1, 1)), m)
})

test_that("baseline argument", {
  m <- "should be a single distr\\(\\) specification or character string naming a study in the network, a list of such specifications, or NULL"
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = 1), m)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = NA), m)

  m2 <- "or a list of length 2"
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(distr(qnorm, 1, 1),
                                                                                  distr(qnorm, 2, 1),
                                                                                  distr(qnorm, 3, 1))), m2)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(One = distr(qnorm, 1, 1),
                                                                                  Two = distr(qnorm, 2, 1),
                                                                                  Three = distr(qnorm, 3, 1))), m2)

  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(One = distr(qnorm, 1, 1),
                                                                                  Three = distr(qnorm, 3, 1))),
               "must match all study names")

  m3 <- "must match the name of an IPD or AgD \\(arm-based\\) study in the network"
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = "a"), m3)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list("a")), m3)
})

test_that("newdata validation", {
  expect_error(predict(pso_fit, newdata = pso_new[-1], study = study, baseline = distr(qnorm, 1, 1)),
               'Regression variable "durnpso" not found in `newdata`')
  expect_error(predict(pso_fit, newdata = pso_new[-(1:2)], study = study, baseline = distr(qnorm, 1, 1)),
               'Regression variables "durnpso" and "prevsys" not found in `newdata`')

  make_bad <- function(x, vars = "durnpso") {
    bad <- pso_new
    bad[vars] <- x
    return(bad)
  }

  expect_error(predict(pso_fit, newdata = make_bad(NA), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(predict(pso_fit, newdata = make_bad(NaN), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(predict(pso_fit, newdata = make_bad(Inf), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")

  expect_error(predict(pso_fit, newdata = make_bad(NA, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(predict(pso_fit, newdata = make_bad(NaN, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(predict(pso_fit, newdata = make_bad(Inf, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
})

test_that(".study, .trt columns are correct", {
  pred1 <- tibble::as_tibble(predict(pso_fit))
  expect_identical(paste0("pred[", pred1$.study, ": ", pred1$.trt, "]"),
                   pred1$parameter)

  pred2 <- tibble::as_tibble(predict(pso_fit, newdata = pso_new, study = study, baseline = distr(qnorm, 1, 1)))
  expect_identical(paste0("pred[", pred2$.study, ": ", pred2$.trt, "]"),
                   pred2$parameter)
})


test_that("expand TRUE and FALSE", {
  pred1 <- tibble::as_tibble(predict(pso_fit, expand = TRUE))
  expect_identical(pred1$parameter,
                   paste0("pred[",
                          tidyr::expand_grid(s = pso_net$studies, t = pso_net$treatments) %>%
                            dplyr::mutate(l = paste(s, t, sep = ": ")) %>%
                            dplyr::pull(l), "]"))

  pred2 <- tibble::as_tibble(predict(pso_fit, expand = FALSE))
  expect_identical(pred2$parameter,
                   paste0("pred[",
                          dplyr::distinct(pso_net$ipd, .study, .trt) %>%
                            dplyr::arrange(.study, .trt) %>%
                            dplyr::mutate(l = paste(.study, .trt, sep = ": ")) %>%
                            dplyr::pull(l), "]"))

  pred3 <- tibble::as_tibble(predict(pso_fit, expand = FALSE, baseline = "UNCOVER-1",
                                     newdata = dplyr::mutate(pso_new, .trt = "IXE_Q2W")))
  expect_identical(pred3$parameter, "pred[New 1: IXE_Q2W]")

  expect_error(predict(pso_fit, expand = FALSE, baseline = distr(qnorm, 0, 1), newdata = pso_new),
               "`newdata` should have a `\\.trt` column")

  pred4 <- tibble::as_tibble(predict(pso_fit, expand = FALSE,
                                     baseline = distr(qnorm, 0, 1),
                                     study = study,
                                     newdata = dplyr::mutate(pso_new, .trt = "IXE_Q2W")))
  expect_identical(pred4$parameter,
                   paste0("pred[", c("One", "Two"), ": ", "IXE_Q2W", "]"))

  expect_error(predict(pso_fit, expand = FALSE,
                       baseline = distr(qnorm, 0, 1),
                       newdata = dplyr::mutate(pso_new, .trt = "a")),
               "Treatments in `newdata` do not match those in the network")

  pred5 <- tibble::as_tibble(predict(pso_fit, expand = TRUE, level = "individual"))
  expect_identical(pred5$parameter,
                   paste0("pred[",
                          dplyr::group_by(pso_net$ipd, .study) %>%
                            dplyr::mutate(id = 1:dplyr::n()) %>%
                            dplyr::select(-.trt) %>%
                            tidyr::expand(id, .trt = pso_net$treatments) %>%
                            dplyr::mutate(l = paste0(.study, ": ", .trt, ", ", id)) %>%
                            dplyr::pull(l), "]"))

  pred6 <- tibble::as_tibble(predict(pso_fit, expand = FALSE, level = "individual"))
  expect_identical(pred6$parameter,
                   paste0("pred[", dplyr::group_by(pso_net$ipd, .study, .trt) %>%
                            dplyr::mutate(l = paste0(.study, ": ", .trt, ", ", 1:dplyr::n())) %>%
                            dplyr::pull(l), "]"))
})


hta_net <- set_agd_arm(hta_psoriasis,
                       study = paste(studyc, year),
                       trt = trtc,
                       r = multi(r0 = sample_size - rowSums(cbind(PASI50, PASI75, PASI90), na.rm = TRUE),
                                 PASI50, PASI75, PASI90,
                                 inclusive = FALSE,
                                 type = "ordered"))

# Only small number of samples to test
hta_fit_FE <- suppressWarnings(nma(hta_net,
                               trt_effects = "fixed",
                               link = "probit",
                               prior_intercept = normal(scale = 100),
                               prior_trt = normal(scale = 10),
                               prior_aux = flat(),
                               iter = 10))

test_that(".study, .trt, .category columns are correct", {
  pred1 <- tibble::as_tibble(predict(hta_fit_FE))
  expect_identical(paste0("pred[", pred1$.study, ": ", pred1$.trt, ", ", pred1$.category, "]"),
                   pred1$parameter)

  pred2 <- tibble::as_tibble(predict(hta_fit_FE, baseline = distr(qnorm, 0, 1)))
  expect_identical(paste0("pred[", pred2$.trt, ", ", pred2$.category, "]"),
                   pred2$parameter)
})

ndmm_net <- combine_network(
  set_ipd(dplyr::slice_sample(ndmm_ipd, n = 10, by = c("study", "trt")),
          study, trt,
          Surv = Surv(eventtime/7, status),
          trt_class = trt != "Pbo"),
  set_agd_surv(dplyr::slice_sample(ndmm_agd, n= 10, by = c("study", "trt")),
               study, trt,
               Surv = Surv(eventtime/7, status),
               covariates = ndmm_agd_covs,
               trt_class = trt != "Pbo")
) %>%
  add_integration(age = distr(qnorm, age_mean, age_sd),
                  n_int = 5)

# Only small number of samples to test
ndmm_fit_weib <- suppressWarnings(nma(ndmm_net,
                                  likelihood = "weibull-aft",
                                  prior_intercept = normal(0, 100),
                                  prior_trt = normal(0, 10),
                                  prior_aux = half_normal(10),
                                  iter = 10))

ndmm_fit_weib_nphr <- suppressWarnings(nma(ndmm_net,
                                      likelihood = "weibull-aft",
                                      aux_regression = ~.trt,
                                      prior_intercept = normal(0, 100),
                                      prior_trt = normal(0, 10),
                                      prior_aux = half_normal(10),
                                      prior_aux_reg = normal(0, 10),
                                      iter = 10))

ndmm_fit_weib_nphs <- suppressWarnings(nma(ndmm_net,
                                           likelihood = "weibull-aft",
                                           aux_by = ".trt",
                                           prior_intercept = normal(0, 100),
                                           prior_trt = normal(0, 10),
                                           prior_aux = half_normal(10),
                                           iter = 10))

ndmm_fit_exp <- suppressWarnings(nma(ndmm_net,
                                     likelihood = "exponential",
                                     prior_intercept = normal(0, 100),
                                     prior_trt = normal(0, 10),
                                     iter = 10))

ndmm_fit_gengamma <- suppressWarnings(nma(ndmm_net,
                                          likelihood = "gengamma",
                                          prior_intercept = normal(0, 100),
                                          prior_trt = normal(0, 10),
                                          prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                          iter = 10))

ndmm_fit_gengamma_nphr <- suppressWarnings(nma(ndmm_net,
                                          likelihood = "gengamma",
                                          aux_regression = ~.trt,
                                          prior_intercept = normal(0, 100),
                                          prior_trt = normal(0, 10),
                                          prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                          prior_aux_reg = normal(0, 10),
                                          iter = 10))

ndmm_fit_gengamma_nphs <- suppressWarnings(nma(ndmm_net,
                                               likelihood = "gengamma",
                                               aux_by = ".trt",
                                               prior_intercept = normal(0, 100),
                                               prior_trt = normal(0, 10),
                                               prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                               iter = 10))

ndmm_fit_mspline <- suppressWarnings(nma(ndmm_net,
                                         likelihood = "mspline",
                                         prior_intercept = normal(0, 100),
                                         prior_trt = normal(0, 10),
                                         prior_aux = half_normal(1),
                                         iter = 10))

ndmm_fit_mspline_nphr <- suppressWarnings(nma(ndmm_net,
                                         likelihood = "mspline",
                                         aux_regression = ~.trt,
                                         prior_intercept = normal(0, 100),
                                         prior_trt = normal(0, 10),
                                         prior_aux = half_normal(1),
                                         prior_aux_reg = normal(0, 10),
                                         iter = 10))

ndmm_fit_mspline_nphs <- suppressWarnings(nma(ndmm_net,
                                              likelihood = "mspline",
                                              aux_by = ".trt",
                                              prior_intercept = normal(0, 100),
                                              prior_trt = normal(0, 10),
                                              prior_aux = half_normal(1),
                                              iter = 10))

ndmm_preddat <- dplyr::bind_rows(
  dplyr::transmute(ndmm_net$ipd, .study, .trt, .time = .Surv[, "time"]),
  tidyr::unnest(ndmm_net$agd_arm, cols = ".Surv") %>% dplyr::transmute(.study, .trt, .time = .Surv[, "time"]),
)

test_that(".study, .trt, .time columns are correct (weibull, no regression, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))

  # Check non-PH
  pred1.1_nphr <- tibble::as_tibble(predict(ndmm_fit_weib_nphr, type = "survival"))
  expect_equivalent(pred1.1_nphr[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphr$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::inner_join(preddat1,
                                     dplyr::distinct(ndmm_preddat, .study, .trt),
                                     by = c(".study", ".trt")) %>%
    dplyr::arrange(.study, .trt)

  pred1.1_nphs <- tibble::as_tibble(predict(ndmm_fit_weib_nphs, type = "survival"))
  expect_equivalent(pred1.1_nphs[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphs$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

})


test_that(".study, .trt, .time columns are correct (exponential, no regression, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (gengamma, no regression, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))

  # Check non-PH
  pred1.1_nphr <- tibble::as_tibble(predict(ndmm_fit_gengamma_nphr, type = "survival"))
  expect_equivalent(pred1.1_nphr[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphr$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::inner_join(preddat1,
                                     dplyr::distinct(ndmm_preddat, .study, .trt),
                                     by = c(".study", ".trt")) %>%
    dplyr::arrange(.study, .trt)

  pred1.1_nphs <- tibble::as_tibble(predict(ndmm_fit_gengamma_nphs, type = "survival"))
  expect_equivalent(pred1.1_nphs[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphs$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

})


test_that(".study, .trt, .time columns are correct (mspline, no regression, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # expect_warning(
  #   pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "mean")),
  #   "Evaluating M-spline at times beyond the boundary knots")
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  expect_warning(
    pred3.2 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "median")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  expect_warning(
    pred4.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "quantile", quantiles = qs)),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))

  # Check non-PH
  pred1.1_nphr <- tibble::as_tibble(predict(ndmm_fit_mspline_nphr, type = "survival"))
  expect_equivalent(pred1.1_nphr[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphr$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::inner_join(preddat1,
                                     dplyr::distinct(ndmm_preddat, .study, .trt),
                                     by = c(".study", ".trt")) %>%
    dplyr::arrange(.study, .trt)

  pred1.1_nphs <- tibble::as_tibble(predict(ndmm_fit_mspline_nphs, type = "survival"))
  expect_equivalent(pred1.1_nphs[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1_nphs$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

})

test_that("Survival predictions for new studies require correct args", {
  m <- "Specify both `baseline` and `aux`, or neither"

  expect_error(predict(ndmm_fit_weib, type = "survival",
                       times = 0:5,
                       baseline = distr(qnorm, 0, 1)),
               m)
  expect_error(predict(ndmm_fit_weib, type = "survival",
                       times = 0:5,
                       aux = distr(qnorm, 0, 1)),
               m)
  expect_error(predict(ndmm_fit_weib, type = "survival",
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qnorm, 0, 1)),
               "`times` must be specified")
})

test_that(".study, .trt, .time columns are correct (weibull, no regression, new data)", {

  time <- 0:5

  # Prediction format new times
  preddat1 <- tidyr::expand_grid(.study = "New 1",
                                 .trt = unique(ndmm_preddat$.trt),
                                 .time = time) %>%
    dplyr::group_by(.trt) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "survival",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.1[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_weib_nphr, type = "survival",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.1r[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1s <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_weib_nphs, type = "survival",
                                        times = time,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".trt", ".time")],
                    preddat1s[, c(".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1s$.trt, ", ", preddat1s$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "hazard",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.2[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "cumhaz",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.3[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tibble::tibble(.trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "mean",
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = distr(qlnorm, 0, 0.01)))
  # expect_equivalent(pred3.1[, ".trt"],
  #                   preddat3)
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "median",
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.2[, ".trt"],
                    preddat3)
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "rmst",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.3[, c(".trt", ".time")],
                    dplyr::mutate(preddat3, .time = 3))
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "link",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.4[, ".trt"],
                    preddat3)
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- tidyr::expand_grid(.trt = unique(ndmm_preddat$.trt),
                                 .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib, type = "quantile", quantiles = qs,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (exponential, no regression, new data)", {

  time <- 0:5

  # Prediction format new times
  preddat1 <- tidyr::expand_grid(.study = "New 1",
                                 .trt = unique(ndmm_preddat$.trt),
                                 .time = time) %>%
    dplyr::group_by(.trt) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "survival",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.1[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "hazard",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.2[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "cumhaz",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.3[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tibble::tibble(.trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "mean",
  #                                      baseline = distr(qnorm, 0, 1)))
  # expect_equivalent(pred3.1[, ".trt"],
  #                   preddat3)
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "median",
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.2[, ".trt"],
                    preddat3)
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "rmst",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.3[, c(".trt", ".time")],
                    dplyr::mutate(preddat3, .time = 3))
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "link",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.4[, ".trt"],
                    preddat3)
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- tidyr::expand_grid(.trt = unique(ndmm_preddat$.trt),
                                 .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp, type = "quantile", quantiles = qs,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (gengamma, no regression, new data)", {

  time <- 0:5

  # Prediction format new times
  preddat1 <- tidyr::expand_grid(.study = "New 1",
                                 .trt = unique(ndmm_preddat$.trt),
                                 .time = time) %>%
    dplyr::group_by(.trt) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "survival",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.1[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_gengamma_nphr, type = "survival",
                                        times = time,
                                        baseline = distr(qnorm, 0, 1),
                                        aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                   k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.1r[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1s <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_gengamma_nphs, type = "survival",
                                        times = time,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".trt", ".time")],
                    preddat1s[, c(".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1s$.trt, ", ", preddat1s$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "hazard",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.2[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "cumhaz",
                                       times = time,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.3[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tibble::tibble(.trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "mean",
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = list(sigma = distr(qlnorm, 0, 0.1),
  #                                                 k = distr(qlnorm, 0, 0.1))))
  # expect_equivalent(pred3.1[, ".trt"],
  #                   preddat3)
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "median",
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.2[, ".trt"],
                    preddat3)
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "rmst",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.3[, c(".trt", ".time")],
                    dplyr::mutate(preddat3, .time = 3))
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "link",
                                       time = 3,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.4[, ".trt"],
                    preddat3)
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- tidyr::expand_grid(.trt = unique(ndmm_preddat$.trt),
                                 .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma, type = "quantile", quantiles = qs,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that("errors for aux lists", {
  expect_error(predict(ndmm_fit_gengamma, type = "survival",
                       times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = list(k = distr(qlnorm, 0, 0.1))),
               "`aux` must be a named list of distr\\(\\) specifications for sigma and k")
  expect_error(predict(ndmm_fit_gengamma, type = "survival",
                       times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 0.1))),
               "`aux` must be a named list of distr\\(\\) specifications for sigma and k")
  expect_error(predict(ndmm_fit_gengamma, type = "survival",
                       times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qlnorm, 0, 0.1)),
               "`aux` must be a named list of distr\\(\\) specifications for sigma and k")
})

test_that(".study, .trt, .time columns are correct (mspline, no regression, new data)", {

  expect_error(predict(ndmm_fit_mspline, type = "survival",
                       times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qnorm, 0, 1)),
               'Producing predictions with external `aux` spline coefficients is not currently supported for "mspline" models.')

  time <- 0:5

  # Prediction format new times
  preddat1 <- tidyr::expand_grid(.study = "New 1",
                                 .trt = unique(ndmm_preddat$.trt),
                                 .time = time) %>%
    dplyr::group_by(.trt) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "survival",
                                       times = time,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred1.1[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_mspline_nphr, type = "survival",
                                        times = time,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1r[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1s <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_mspline_nphs, type = "survival",
                                        times = time,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".trt", ".time")],
                    preddat1s[, c(".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1s$.trt, ", ", preddat1s$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "hazard",
                                       times = time,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred1.2[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "cumhaz",
                                       times = time,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred1.3[, c(".trt", ".time")],
                    preddat1[, c(".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tibble::tibble(.trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "mean",
  #                                      baseline = "Attal2012,
  #                                      aux = "Attal2012"))
  # expect_equivalent(pred3.1[, ".trt"],
  #                   preddat3)
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.trt, "]"))

  pred3.2 <- suppressWarnings(tibble::as_tibble(predict(ndmm_fit_mspline, type = "median",
                                       baseline = "Attal2012",
                                       aux = "Attal2012")))
  expect_equivalent(pred3.2[, ".trt"],
                    preddat3)
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "rmst",
                                       time = 3,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred3.3[, c(".trt", ".time")],
                    dplyr::mutate(preddat3, .time = 3))
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline, type = "link",
                                       time = 3,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred3.4[, ".trt"],
                    preddat3)
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- tidyr::expand_grid(.trt = unique(ndmm_preddat$.trt),
                                 .quantile = qs)

  pred4.1 <- suppressWarnings(tibble::as_tibble(predict(ndmm_fit_mspline, type = "quantile", quantiles = qs,
                                       baseline = "Attal2012",
                                       aux = "Attal2012")))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})


ndmm_fit_weib_reg <- suppressWarnings(nma(ndmm_net,
                                          likelihood = "weibull-aft",
                                          regression = ~I(age/10)*.trt,
                                          class_interactions= "common",
                                          prior_intercept = normal(0, 100),
                                          prior_trt = normal(0, 10),
                                          prior_reg = normal(0, 10),
                                          prior_aux = half_normal(10),
                                          iter = 10,
                                          seed = 42))

ndmm_fit_weib_reg_nphr <- suppressWarnings(nma(ndmm_net,
                                          likelihood = "weibull-aft",
                                          regression = ~I(age/10)*.trt,
                                          aux_regression = ~.trt,
                                          class_interactions= "common",
                                          prior_intercept = normal(0, 100),
                                          prior_trt = normal(0, 10),
                                          prior_reg = normal(0, 10),
                                          prior_aux = half_normal(10),
                                          iter = 10,
                                          seed = 42))

ndmm_fit_weib_reg_nphs <- suppressWarnings(nma(ndmm_net,
                                          likelihood = "weibull-aft",
                                          regression = ~I(age/10)*.trt,
                                          aux_by = ".trt",
                                          class_interactions= "common",
                                          prior_intercept = normal(0, 100),
                                          prior_trt = normal(0, 10),
                                          prior_reg = normal(0, 10),
                                          prior_aux = half_normal(10),
                                          iter = 10,
                                          seed = 42))

ndmm_fit_exp_reg <- suppressWarnings(nma(ndmm_net,
                                         likelihood = "exponential",
                                         regression = ~I(age/10)*.trt,
                                         class_interactions = "common",
                                         prior_intercept = normal(0, 100),
                                         prior_trt = normal(0, 10),
                                         prior_reg = normal(0, 10),
                                         iter = 10))

ndmm_fit_gengamma_reg <- suppressWarnings(nma(ndmm_net,
                                              likelihood = "gengamma",
                                              regression = ~I(age/10)*.trt,
                                              class_interactions= "common",
                                              prior_intercept = normal(0, 100),
                                              prior_trt = normal(0, 10),
                                              prior_reg = normal(0, 10),
                                              prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                              init_r = 0.1,
                                              iter = 10))

ndmm_fit_gengamma_reg_nphr <- suppressWarnings(nma(ndmm_net,
                                              likelihood = "gengamma",
                                              regression = ~I(age/10)*.trt,
                                              aux_regression = ~.trt,
                                              class_interactions= "common",
                                              prior_intercept = normal(0, 100),
                                              prior_trt = normal(0, 10),
                                              prior_reg = normal(0, 10),
                                              prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                              init_r = 0.1,
                                              iter = 10))

ndmm_fit_gengamma_reg_nphs <- suppressWarnings(nma(ndmm_net,
                                              likelihood = "gengamma",
                                              regression = ~I(age/10)*.trt,
                                              aux_by = ".trt",
                                              class_interactions= "common",
                                              prior_intercept = normal(0, 100),
                                              prior_trt = normal(0, 10),
                                              prior_reg = normal(0, 10),
                                              prior_aux = list(sigma = half_normal(5), k = half_normal(5)),
                                              init_r = 0.1,
                                              iter = 10))

ndmm_fit_mspline_reg <- suppressWarnings(nma(ndmm_net,
                                             likelihood = "mspline",
                                             regression = ~I(age/10)*.trt,
                                             class_interactions = "common",
                                             prior_intercept = normal(0, 100),
                                             prior_trt = normal(0, 10),
                                             prior_reg = normal(0, 10),
                                             prior_aux = half_normal(1),
                                             iter = 10))

ndmm_fit_mspline_reg_nphr <- suppressWarnings(nma(ndmm_net,
                                             likelihood = "mspline",
                                             regression = ~I(age/10)*.trt,
                                             aux_regression = ~.trt,
                                             class_interactions = "common",
                                             prior_intercept = normal(0, 100),
                                             prior_trt = normal(0, 10),
                                             prior_reg = normal(0, 10),
                                             prior_aux = half_normal(1),
                                             iter = 10))

ndmm_fit_mspline_reg_nphs <- suppressWarnings(nma(ndmm_net,
                                             likelihood = "mspline",
                                             regression = ~I(age/10)*.trt,
                                             aux_by = ".trt",
                                             class_interactions = "common",
                                             prior_intercept = normal(0, 100),
                                             prior_trt = normal(0, 10),
                                             prior_reg = normal(0, 10),
                                             prior_aux = half_normal(1),
                                             iter = 10))


test_that("Survival predictions for new studies require correct args", {
  m <- "Specify all of `newdata`, `baseline`, and `aux`, or none"

  expect_error(predict(ndmm_fit_weib_reg, type = "survival",
                       baseline = distr(qnorm, 0, 1)),
               m)
  expect_error(predict(ndmm_fit_weib_reg, type = "survival",
                       times = 0:5,
                       aux = distr(qnorm, 0, 1)),
               m)
  expect_error(predict(ndmm_fit_weib_reg, type = "survival",
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qnorm, 0, 1)),
               m)
})

test_that(".study, .trt, .time columns are correct (weibull, regression, individual, network data)", {
  # Prediction format for level = "individual", observed times
  preddat1 <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt))

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", level = "individual"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphr, type = "survival", level = "individual"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphs, type = "survival", level = "individual"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", level = "individual"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", level = "individual"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Single summaries, also at individual times
  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "mean", level = "individual"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "median", level = "individual"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst", level = "individual"))
  # expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.3$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "link", level = "individual"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat1, tibble::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "quantile", quantiles = qs,
                                       level = "individual"))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (exponential, regression, individual, network data)", {
  # Prediction format for level = "individual", observed times
  preddat1 <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt))

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", level = "individual"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", level = "individual"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", level = "individual"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Single summaries, also at individual times
  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "mean", level = "individual"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "median", level = "individual"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst", level = "individual"))
  # expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.3$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "link", level = "individual"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat1, tibble::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "quantile", quantiles = qs,
                                       level = "individual"))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (gengamma, regression, individual, network data)", {
  # Prediction format for level = "individual", observed times
  preddat1 <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt))

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", level = "individual"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphr, type = "survival", level = "individual"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphs, type = "survival", level = "individual"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", level = "individual"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", level = "individual"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Single summaries, also at individual times
  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "mean", level = "individual"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "median", level = "individual"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst", level = "individual"))
  # expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.3$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "link", level = "individual"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat1, tibble::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "quantile", quantiles = qs,
                                       level = "individual"))
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (mspline, regression, individual, network data)", {
  # Prediction format for level = "individual", observed times
  preddat1 <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt))

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", level = "individual"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphr, type = "survival", level = "individual"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::filter(ndmm_preddat, .study %in% unique(ndmm_ipd$study)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup()

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphs, type = "survival", level = "individual"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", level = "individual"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", level = "individual"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Single summaries, also at individual times
  # expect_warning(
  #   pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "mean", level = "individual")),
  #   "Evaluating M-spline at times beyond the boundary knots")
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  expect_warning(
    pred3.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "median", level = "individual")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst", level = "individual"))
  # expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                   preddat1[, c(".study", ".trt")])
  # expect_identical(pred3.3$parameter,
  #                  paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "link", level = "individual"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat1[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat1, tibble::tibble(.quantile = qs))

  expect_warning(
    pred4.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "quantile", quantiles = qs,
                                        level = "individual")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred4.1[, c(".trt", ".quantile")],
                    preddat4[, c(".trt", ".quantile")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (weibull, regression, aggregate, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphr, type = "survival"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.study, .trt, id)

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphs, type = "survival"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (exponential, regression, aggregate, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (gengamma, regression, aggregate, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphr, type = "survival"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.study, .trt, id)

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphs, type = "survival"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "median"))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "quantile", quantiles = qs))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (mspline, regression, aggregate, network data)", {
  # Prediction format for observed times
  preddat1 <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.trt) %>%
    dplyr::cross_join(dplyr::distinct(ndmm_preddat, .trt)) %>%
    dplyr::arrange(.study, .trt)

  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphr, type = "survival"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- ndmm_preddat %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.study, .trt, id)

  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphs, type = "survival"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for provided times
  times <- 0:5
  preddat2 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt) %>%
    dplyr::mutate(.time = list(times), id = list(seq_along(times))) %>%
    tidyr::unnest(cols = c(".time", "id"))

  pred2.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", times = times))
  expect_equivalent(pred2.1[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.1$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", times = times))
  expect_equivalent(pred2.2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", times = times))
  expect_equivalent(pred2.3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt)

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "mean"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  expect_warning(
    pred3.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "median")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "rmst"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "link"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::distinct(ndmm_preddat, .study, .trt) %>%
    tidyr::expand(.study, .trt, .quantile = qs)

  expect_warning(
    pred4.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "quantile", quantiles = qs)),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (weibull, regression, aggregate, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm) %>%
    add_integration(age = distr(qlnorm, meanlog = log(4.5), sdlog = 1), n_int = 5)

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphr, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  newdata2 <- dplyr::bind_rows(
    dplyr::mutate(newdata, study = "Test"),
    dplyr::mutate(newdata, study = "B"),
    dplyr::mutate(newdata, study = "C")
  )
  preddat2 <- dplyr::mutate(newdata2, .study = factor(study), .time = rep(tm, times = 3)) %>%
    dplyr::group_by(.study) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  pred1.1r2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphr, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata2,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1r2[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r2$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred1.1r3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphr, type = "survival", time = time,
                                         study = study,
                                         newdata = newdata2,
                                         baseline = list(Test = "Attal2012", B = "Morgan2012", C = distr(qlnorm, 0, 1)),
                                         aux = list(Test = "Attal2012", B = "Morgan2012", C = distr(qlnorm, 0, 0.01))))
  expect_equivalent(pred1.1r3[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r3$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  preddat1_nphs <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphs, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  preddat2_nphs <- dplyr::filter(preddat2, .trt %in% c("Pbo", "Len"))
  pred1.1s2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphs, type = "survival", time = time,
                                         study = study,
                                         newdata = newdata2,
                                         baseline = "Attal2012",
                                         aux = "Attal2012"))
  expect_equivalent(pred1.1s2[, c(".study", ".trt", ".time")],
                    preddat2_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s2$parameter,
                   paste0("pred[", preddat2_nphs$.study, ": ", preddat2_nphs$.trt, ", ", preddat2_nphs$id, "]"))

  preddat3_nphs <- dplyr::filter(preddat2, .study == "Test" & .trt %in% c("Pbo", "Len") | .study == "B" & .trt %in% c("Pbo", "Thal") | .study == "C")
  pred1.1s3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg_nphs, type = "survival", time = time,
                                         study = study,
                                         newdata = newdata2,
                                         baseline = list(Test = "Attal2012", B = "Morgan2012", C = distr(qlnorm, 0, 1)),
                                         aux = list(Test = "Attal2012", B = "Morgan2012", C = distr(qlnorm, 0, 0.01))))
  expect_equivalent(pred1.1s3[, c(".study", ".trt", ".time")],
                    preddat3_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s3$parameter,
                   paste0("pred[", preddat3_nphs$.study, ": ", preddat3_nphs$.trt, ", ", preddat3_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  pred1.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.4[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.5 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.5[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.5$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.6 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.6[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.6$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tidyr::expand_grid(.study = unique(factor(newdata$study)),
                                 .trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = distr(qlnorm, 0, 0.01)))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # pred3.3 <- withr::with_seed(
  #   123,
  #   tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "rmst", time = 3,
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = distr(qlnorm, 0, 0.01)))
  # )

  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- preddat3 %>% tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "quantile", quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})


test_that(".study, .trt, .time columns are correct (exponential, regression, aggregate, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm) %>%
    add_integration(age = distr(qlnorm, meanlog = log(4.5), sdlog = 1), n_int = 5)

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  pred1.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.4[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.5 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.5[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.5$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.6 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.6[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.6$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tidyr::expand_grid(.study = unique(factor(newdata$study)),
                                 .trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1)))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "rmst", time = 3,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- preddat3 %>% tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "quantile", quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (gengamma, regression, aggregate, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm) %>%
    add_integration(age = distr(qlnorm, meanlog = log(4.5), sdlog = 1), n_int = 5)

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphr, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg_nphs, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))


  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  pred1.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.4[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.5 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.5[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.5$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.6 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.6[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.6$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tidyr::expand_grid(.study = unique(factor(newdata$study)),
                                 .trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = list(sigma = distr(qlnorm, 0, 0.1),
  #                                                 k = distr(qlnorm, 0, 0.1))))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  #pred3.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "rmst", time = 3,
  #                                     study = study,
  #                                     newdata = newdata,
  #                                     baseline = distr(qnorm, 0, 1),
  #                                     aux = list(sigma = distr(qlnorm, 0, 0.1),
  #                                                k = distr(qlnorm, 0, 0.1))))
  #expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                  preddat3[, c(".study", ".trt")])
  #expect_identical(pred3.3$parameter,
  #                 paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- preddat3 %>% tidyr::expand(.study, .trt, .quantile = qs)

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "quantile", quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (mspline, regression, aggregate, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm) %>%
    add_integration(age = distr(qlnorm, meanlog = log(4.5), sdlog = 1), n_int = 5)

  expect_error(predict(ndmm_fit_mspline_reg, type = "survival", time = time,
                       study = study,
                       newdata = newdata,
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qlnorm, 0, 0.01)),
               'Producing predictions with external `aux` spline coefficients is not currently supported for "mspline" models.')

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1b <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred1.1b[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1b$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1r <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphr, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1r[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1r$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  preddat1_nphs <- dplyr::filter(preddat1, .trt %in% c("Pbo", "Len"))
  pred1.1s <- tibble::as_tibble(predict(ndmm_fit_mspline_reg_nphs, type = "survival", time = time,
                                        study = study,
                                        newdata = newdata,
                                        baseline = "Attal2012",
                                        aux = "Attal2012"))
  expect_equivalent(pred1.1s[, c(".study", ".trt", ".time")],
                    preddat1_nphs[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1s$parameter,
                   paste0("pred[", preddat1_nphs$.study, ": ", preddat1_nphs$.trt, ", ", preddat1_nphs$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  pred1.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.4[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.4$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.5 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.5[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.5$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.6 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.6[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.6$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Prediction format for single summaries
  preddat3 <- tidyr::expand_grid(.study = unique(factor(newdata$study)),
                                 .trt = unique(ndmm_preddat$.trt))

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = "Attal2012"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  expect_warning(
    pred3.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "median",
                                         study = study,
                                         newdata = newdata,
                                         baseline = distr(qnorm, 0, 1),
                                         aux = "Attal2012")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # pred3.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "rmst", time = 3,
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = "Attal2012"))
  # expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.3$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- preddat3 %>% tidyr::expand(.study, .trt, .quantile = qs)

  expect_warning(
    pred4.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "quantile", quantiles = qs,
                                         study = study,
                                         newdata = newdata,
                                         baseline = distr(qnorm, 0, 1),
                                         aux = "Attal2012")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (weibull, regression, individual, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm, age = rlnorm(length(tm), log(4.5), 0.25))

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  # With external times vector, evaluate every time and every treatment for every new individual
  preddat2 <- dplyr::mutate(newdata, .study = factor(study)) %>%
    dplyr::cross_join(
      dplyr::cross_join(dplyr::tibble(.time = tm, id = 1:length(tm)),
                        dplyr::tibble(.trt = unique(ndmm_preddat$.trt))))
  pred2.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "survival", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred2.4[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.4$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.5 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "hazard", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred2.5[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.5$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.6 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "cumhaz", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred2.6[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.6$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- preddat1

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, level = "individual",
  #                                      type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = distr(qlnorm, 0, 0.01)))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, level = "individual",
                                       type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, level = "individual",
                                       type = "rmst", time = 3,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, level = "individual",
                                       type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat3, dplyr::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_weib_reg, type = "quantile", level = "individual",
                                       quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = distr(qlnorm, 0, 0.01)))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (exponential, regression, individual, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm, age = rlnorm(length(tm), log(4.5), 0.25))

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  # With external times vector, evaluate every time and every treatment for every new individual
  preddat2 <- dplyr::mutate(newdata, .study = factor(study)) %>%
    dplyr::cross_join(
      dplyr::cross_join(dplyr::tibble(.time = tm, id = 1:length(tm)),
                        dplyr::tibble(.trt = unique(ndmm_preddat$.trt))))
  pred2.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "survival", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred2.4[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.4$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.5 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "hazard", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred2.5[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.5$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.6 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "cumhaz", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred2.6[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.6$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- preddat1

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, level = "individual",
  #                                      type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1)))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, level = "individual",
                                       type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, level = "individual",
                                       type = "rmst", time = 3,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, level = "individual",
                                       type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat3, dplyr::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_exp_reg, type = "quantile", level = "individual",
                                       quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1)))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (gengamma, regression, individual, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm, age = rlnorm(length(tm), log(4.5), 0.25))

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  # With external times vector, evaluate every time and every treatment for every new individual
  preddat2 <- dplyr::mutate(newdata, .study = factor(study)) %>%
    dplyr::cross_join(
      dplyr::cross_join(dplyr::tibble(.time = tm, id = 1:length(tm)),
                        dplyr::tibble(.trt = unique(ndmm_preddat$.trt))))
  pred2.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "survival", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred2.4[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.4$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.5 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "hazard", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred2.5[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.5$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.6 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "cumhaz", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred2.6[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.6$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- preddat1

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, level = "individual",
  #                                      type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = list(sigma = distr(qlnorm, 0, 0.1),
  #                                                 k = distr(qlnorm, 0, 0.1))))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.2 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, level = "individual",
                                       type = "median",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  #pred3.3 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, level = "individual",
  #                                     type = "rmst", time = 3,
  #                                     study = study,
  #                                     newdata = newdata,
  #                                     baseline = distr(qnorm, 0, 1),
  #                                     aux = list(sigma = distr(qlnorm, 0, 0.1),
  #                                                k = distr(qlnorm, 0, 0.1))))
  #expect_equivalent(pred3.3[, c(".study", ".trt")],
  #                  preddat3[, c(".study", ".trt")])
  #expect_identical(pred3.3$parameter,
  #                 paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, level = "individual",
                                       type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat3, dplyr::tibble(.quantile = qs))

  pred4.1 <- tibble::as_tibble(predict(ndmm_fit_gengamma_reg, type = "quantile", level = "individual",
                                       quantiles = qs,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = list(sigma = distr(qlnorm, 0, 0.1),
                                                  k = distr(qlnorm, 0, 0.1))))
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})

test_that(".study, .trt, .time columns are correct (mspline, regression, individual, new data)", {

  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm, age = rlnorm(length(tm), log(4.5), 0.25))

  expect_error(predict(ndmm_fit_mspline_reg, type = "survival", time = time, level = "individual",
                       study = study,
                       newdata = newdata,
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qlnorm, 0, 0.01)),
               'Producing predictions with external `aux` spline coefficients is not currently supported for "mspline" models.')

  # Prediction format new times
  preddat1 <- dplyr::mutate(newdata, .study = factor(study), .time = tm, id = 1:dplyr::n()) %>%
    dplyr::cross_join(dplyr::tibble(.trt = unique(ndmm_preddat$.trt)))

  # Times from newdata column
  pred1.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.1[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.1b <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = "Attal2012",
                                       aux = "Attal2012"))
  expect_equivalent(pred1.1b[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.1b$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.2[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.2$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  pred1.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", level = "individual",
                                       time = time,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred1.3[, c(".study", ".trt", ".time")],
                    preddat1[, c(".study", ".trt", ".time")])
  expect_identical(pred1.3$parameter,
                   paste0("pred[", preddat1$.study, ": ", preddat1$.trt, ", ", preddat1$id, "]"))

  # Times from global env
  # With external times vector, evaluate every time and every treatment for every new individual
  preddat2 <- dplyr::mutate(newdata, .study = factor(study)) %>%
    dplyr::cross_join(
      dplyr::cross_join(dplyr::tibble(.time = tm, id = 1:length(tm)),
                        dplyr::tibble(.trt = unique(ndmm_preddat$.trt))))
  pred2.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "survival", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred2.4[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.4$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.5 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "hazard", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred2.5[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.5$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  pred2.6 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "cumhaz", level = "individual",
                                       time = tm,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred2.6[, c(".study", ".trt", ".time")],
                    preddat2[, c(".study", ".trt", ".time")])
  expect_identical(pred2.6$parameter,
                   paste0("pred[", preddat2$.study, ": ", preddat2$.trt, ", ", preddat2$id, "]"))

  # Prediction format for single summaries
  preddat3 <- preddat1

  # pred3.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, level = "individual",
  #                                      type = "mean",
  #                                      study = study,
  #                                      newdata = newdata,
  #                                      baseline = distr(qnorm, 0, 1),
  #                                      aux = "Attal2012"))
  # expect_equivalent(pred3.1[, c(".study", ".trt")],
  #                   preddat3[, c(".study", ".trt")])
  # expect_identical(pred3.1$parameter,
  #                  paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  expect_warning(
    pred3.2 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, level = "individual",
                                         type = "median",
                                         study = study,
                                         newdata = newdata,
                                         baseline = distr(qnorm, 0, 1),
                                         aux = "Attal2012")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred3.2[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.2$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.3 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, level = "individual",
                                       type = "rmst", time = 3,
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred3.3[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.3$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  pred3.4 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, level = "individual",
                                       type = "link",
                                       study = study,
                                       newdata = newdata,
                                       baseline = distr(qnorm, 0, 1),
                                       aux = "Attal2012"))
  expect_equivalent(pred3.4[, c(".study", ".trt")],
                    preddat3[, c(".study", ".trt")])
  expect_identical(pred3.4$parameter,
                   paste0("pred[", preddat3$.study, ": ", preddat3$.trt, ", ", preddat3$id, "]"))

  # Prediction format for quantiles
  qs <- c(0.2, 0.4, 0.6, 0.8)
  preddat4 <- dplyr::cross_join(preddat3, dplyr::tibble(.quantile = qs))

  expect_warning(
    pred4.1 <- tibble::as_tibble(predict(ndmm_fit_mspline_reg, type = "quantile", level = "individual",
                                         quantiles = qs,
                                         study = study,
                                         newdata = newdata,
                                         baseline = distr(qnorm, 0, 1),
                                         aux = "Attal2012")),
    "Evaluating M-spline at times beyond the boundary knots")
  expect_equivalent(pred4.1[, c(".study", ".trt")],
                    preddat4[, c(".study", ".trt")])
  expect_identical(pred4.1$parameter,
                   paste0("pred[", preddat4$.study, ": ", preddat4$.trt, ", ", preddat4$id, ", ", preddat4$.quantile, "]"))
})

test_that("errors for aux lists (regression, new data)", {
  tm <- 0:5

  newdata <- dplyr::tibble(study = "Test", time = tm) %>%
    add_integration(age = distr(qlnorm, meanlog = log(4.5), sdlog = 1), n_int = 5)

  expect_error(predict(ndmm_fit_gengamma_reg, type = "survival",
                       newdata = newdata, times = time, study = study,
                       baseline = distr(qnorm, 0, 1),
                       aux = list(k = distr(qlnorm, 0, 0.1))),
               "`aux` must be a single named list of distr\\(\\) specifications for sigma and k")
  expect_error(predict(ndmm_fit_gengamma_reg, type = "survival",
                       newdata = newdata, times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 0.1))),
               "`aux` must be a single named list of distr\\(\\) specifications for sigma and k")
  expect_error(predict(ndmm_fit_gengamma_reg, type = "survival",
                       newdata = newdata, times = time,
                       baseline = distr(qnorm, 0, 1),
                       aux = distr(qlnorm, 0, 0.1)),
               "`aux` must be a single named list of distr\\(\\) specifications for sigma and k")
})

test_that("aux argument", {
  newdata <- tibble::tibble(study = c("a", "b"), age = 4, time = 1)

  # Single aux parameter, no regression
  m <- "`aux` must be specified using distr\\(\\), or the name of an IPD or AgD \\(arm-based\\) study in the network"
  expect_error(predict(ndmm_fit_weib, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = 1), m)
  expect_error(predict(ndmm_fit_weib, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = list("a")), m)
  expect_error(predict(ndmm_fit_weib, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = NA), m)

  expect_error(predict(ndmm_fit_weib, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = "a"),
               "`aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a distr\\(\\) distribution")

  # Single aux parameter, regression
  m2 <- "`aux` must be a single distr\\(\\) specification or study name, or a list of length 2 \\(number of `newdata` studies\\)"
  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = 1), m2)
  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = list("a")), m2)
  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = NA), m2)

  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = "a"),
               "All elements of `aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a distr\\(\\) distribution")

  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(a = distr(qnorm, 1, 1),
                                  c = distr(qnorm, 1, 1))),
               "must match all study names from `newdata`")
  expect_error(predict(ndmm_fit_weib_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(a = distr(qnorm, 1, 1),
                                  b = "bad")),
               "All elements of `aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a distr\\(\\) distribution")

  # Multiple aux parameters, no regression
  m3 <- "`aux` must be a named list of distr\\(\\) specifications for sigma and k, or a study name"

  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = 1), m3)
  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = list("a")), m3)
  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = NA), m3)

  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = "a"),
               "`aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a named list of distr\\(\\) specifications for sigma and k")

  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 1))), m3)
  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(k = distr(qlnorm, 0, 1))), m3)
  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(k = distr(qlnorm, 0, 1),
                                  sigma = "bad")), m3)
  expect_error(predict(ndmm_fit_gengamma, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 1),
                                  k = distr(qlnorm, 0, 1),
                                  bad = distr(qlnorm, 0, 1))), m3)

  # Multiple aux parameters, regression
  m4 <- "`aux` must be a single named list of distr\\(\\) specifications for sigma and k, a study name, or a list of length 2 \\(number of `newdata` studies\\) of such lists"
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = 1), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = list("a")), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = NA), m4)

  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1), aux = "a"),
               "All elements of `aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a list of distr\\(\\) distributions")

  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 1))), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(k = distr(qlnorm, 0, 1))), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 1),
                                  k = distr(qlnorm, 0, 1),
                                  bad = distr(qlnorm, 0, 1))), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(
                                  a = list(sigma = distr(qlnorm, 0, 1),
                                           k = distr(qlnorm, 0, 1))
                                 )), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(
                                  a = list(sigma = distr(qlnorm, 0, 1),
                                           k = distr(qlnorm, 0, 1)),
                                  b = list(sigma = distr(qlnorm, 0, 1))
                                 )), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(
                                  a = list(sigma = distr(qlnorm, 0, 1),
                                           k = distr(qlnorm, 0, 1)),
                                  b = list(sigma = distr(qlnorm, 0, 1),
                                           k = distr(qlnorm, 0, 1)),
                                  c = list(sigma = distr(qlnorm, 0, 1),
                                           k = distr(qlnorm, 0, 1))
                       )), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(sigma = distr(qlnorm, 0, 1),
                                  k = "bad")), m4)
  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(
                         a = list(sigma = distr(qlnorm, 0, 1),
                                  k = distr(qlnorm, 0, 1)),
                         b = list(sigma = distr(qlnorm, 0, 1),
                                  k = "bad")
                       )), m4)

  expect_error(predict(ndmm_fit_gengamma_reg, times = 1, study = study, newdata = newdata, baseline = distr(qnorm, 0, 1),
                       aux = list(
                         a = list(sigma = distr(qlnorm, 0, 1),
                                  k = distr(qlnorm, 0, 1)),
                         b = "bad")
                       ), "All elements of `aux` must match the name of an IPD or AgD \\(arm-based\\) study in the network, or be a list of distr\\(\\) distributions")

})

test_that("baseline and aux are escaped correctly", {
  ndmm_net <- set_ipd(ndmm_ipd,
                      study = gsub("(.+)([0-9]{4})", "\\1 (\\2)$", study),
                      trt = trt,
                      Surv = Surv(eventtime/7, status))

  ndmm_fit_gg <- suppressWarnings(nma(ndmm_net,
                  likelihood = "gengamma",
                  aux_regression = ~.trt,
                  init_r = 0.1,
                  iter = 10))

  par <- c("pred[Len, 1]", "pred[Pbo, 1]")

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit_gg, times = 10, type = "survival",
            aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
    par)

  ndmm_fit_weib <- suppressWarnings(nma(ndmm_net,
                                      likelihood = "weibull",
                                      aux_regression = ~.trt,
                                      iter = 10))
  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit_weib, times = 10, type = "survival",
                             aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
    par)

  ndmm_fit_ms <- suppressWarnings(nma(ndmm_net,
                                        likelihood = "mspline",
                                        aux_regression = ~.trt,
                                        iter = 10))
  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit_ms, times = 10, type = "survival",
                             aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
    par)

  par <- c("pred[New 1: Len, 1]", "pred[New 1: Pbo, 1]")
  par2 <- c("pred[(a): Len, 1]", "pred[(a): Pbo, 1]", "pred[(b): Len, 1]", "pred[(b): Pbo, 1]")

  ndmm_fit2_gg <- suppressWarnings(nma(ndmm_net,
                   likelihood = "gengamma",
                   regression = ~I(age/10)*.trt,
                   aux_regression = ~.trt,
                   init_r = 0.1,
                   iter = 10))

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_gg, times = 10, type = "survival",
                             newdata = data.frame(age = 5),
                             aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
    par)

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_gg, times = 10, type = "survival",
                             newdata = data.frame(age = 5, s = c("(a)", "(b)")), study = s,
                             aux = list("(a)" = "Attal (2012)$",
                                        "(b)" = "McCarthy (2012)$"),
                             baseline = list("(a)" = "Attal (2012)$",
                                             "(b)" = "McCarthy (2012)$")))$parameter,
    par2)

  ndmm_fit2_weib <- suppressWarnings(nma(ndmm_net,
                                       likelihood = "weibull",
                                       regression = ~I(age/10)*.trt,
                                       aux_regression = ~.trt,
                                       iter = 10))

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_weib, times = 10, type = "survival",
                             newdata = data.frame(age = 5),
                             aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
    par)

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_weib, times = 10, type = "survival",
                             newdata = data.frame(age = 5, s = c("(a)", "(b)")), study = s,
                             aux = list("(a)" = "Attal (2012)$",
                                        "(b)" = "McCarthy (2012)$"),
                             baseline = list("(a)" = "Attal (2012)$",
                                             "(b)" = "McCarthy (2012)$")))$parameter,
    par2)

  ndmm_fit2_ms <- suppressWarnings(nma(ndmm_net,
                                       likelihood = "mspline",
                                       regression = ~I(age/10)*.trt,
                                       aux_regression = ~.trt,
                                       iter = 10))

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_ms, times = 10, type = "survival",
                             newdata = data.frame(age = 5),
                             aux = "Attal (2012)$", baseline = "Attal (2012)$"))$parameter,
  par)

  expect_identical(
    dplyr::as_tibble(predict(ndmm_fit2_ms, times = 10, type = "survival",
                             newdata = data.frame(age = 5, s = c("(a)", "(b)")), study = s,
                             aux = list("(a)" = "Attal (2012)$",
                                        "(b)" = "McCarthy (2012)$"),
                             baseline = list("(a)" = "Attal (2012)$",
                                             "(b)" = "McCarthy (2012)$")))$parameter,
    par2)

})
