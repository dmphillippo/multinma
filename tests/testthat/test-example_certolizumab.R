test_that("TSD3 Certolizumab example can be reproduced", {
  net <- set_agd_arm(
    certolizumab,
    study = study, trt = trt, n = n, r = r,
    trt_class = dplyr::if_else(trt == "Placebo", "Placebo", "Treatment")
  )

  fit_fe <- nma(
    net,
    regression = ~.mu:.trt,
    prior_intercept = normal(scale = sqrt(1000)),
    prior_trt = normal(scale = 100),
    prior_reg = normal(scale = 100),
    seed = 641
  )

  expect_equal(fit_fe$xbar[[".mu"]], -2.421, tolerance = 0.01)

  summary_fe <- as.data.frame(summary(fit_fe, pars = "beta"))

  expect_equal(summary_fe$mean, -0.93, tolerance = 0.02)
  expect_equal(summary_fe$sd, 0.09, tolerance = 0.02)
  expect_equal(summary_fe$`2.5%`, -1.03, tolerance = 0.06)
  expect_equal(summary_fe$`97.5%`, -0.69, tolerance = 0.06)

  # Note that TSD3 places a uniform [0, 5] prior on the between-trail standard deviation,
  # which is not possible in multinma
  fit_re <- nma(
    net,
    trt_effects = "random",
    regression = ~.mu:.trt,
    prior_intercept = normal(scale = sqrt(1000)),
    prior_trt = normal(scale = 100),
    prior_reg = normal(scale = 100),
    seed = 970
  )

  expect_equal(fit_re$xbar[[".mu"]], -2.421, tolerance = 0.01)

  summary_re <- as.data.frame(summary(fit_re, pars = "beta"))

  expect_equal(summary_re$mean, -0.95, tolerance = 0.02)
  expect_equal(summary_re$sd, 0.10, tolerance = 0.02)
  expect_equal(summary_re$`2.5%`, -1.10, tolerance = 0.06)
  expect_equal(summary_re$`97.5%`, -0.70, tolerance = 0.06)
})
