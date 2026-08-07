
dat <- data.frame(study = 1, trt = c("A", "B"), r = 1, n = 1, x1 = 1)
net_i <- set_ipd(dat, study, trt, r = r)
fit_i <- nma(net_i, regression = ~x1*.trt,
           prior_intercept = normal(0, 100),
           prior_trt = normal(0, 10),
           prior_reg = normal(0, 10),
           test_grad = TRUE)

net_a <- set_agd_arm(dat, study, trt, r = r, n = n)
fit_a <- nma(net_a, regression = ~x1*.trt,
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10),
             test_grad = TRUE)

net_s <- set_ipd(dat, study, trt, Surv = Surv(r, n))
fit_s <- nma(net_s, regression = ~x1*.trt, likelihood = "exponential",
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10),
             test_grad = TRUE)

test_that("argument checks", {
  expect_error(bayes_R2(fit_i, summary = "a"), "must be TRUE or FALSE")
  expect_error(loo_R2(fit_i, summary = "a"), "must be TRUE or FALSE")

  expect_error(bayes_R2(fit_a), "No IPD")
  expect_error(loo_R2(fit_a), "No IPD")

  expect_error(bayes_R2(fit_s), "Not supported for survival outcomes.")
  expect_error(loo_R2(fit_s), "Not supported for survival outcomes.")
})
