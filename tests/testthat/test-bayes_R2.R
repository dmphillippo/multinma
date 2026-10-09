
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


skip_on_cran()

test_that("argument checks", {
  expect_error(loo::loo_predictive_metric(fit_a), "No IPD")
  expect_error(loo::loo_predictive_metric(fit_s), "Not supported for survival outcomes.")
})


# Generate some IPD
set.seed(24601)
n <- 500
dat <- dplyr::tibble(
  study = 1,
  trt = rep(c("A", "B"), each = n/2),
  mu = rnorm(n, 0, 0.2),
  x1 = rnorm(n, 1, 0.5),
  d = ifelse(trt == "B", -1, 0),
  lp = mu + x1 * 0.2 + (trt == "B")*x1*-0.3 + d,
  # continuous
  y = rnorm(n, lp, 0.05),
  # poisson
  E = runif(n, 0.5, 10),
  rp = rpois(n, exp(y) * E),
  # binomial
  r = rbinom(n, 1, plogis(lp)),
  # ordered
  cat = findInterval(y, c(-Inf, 0, 0.5, Inf)),
  r0 = as.numeric(cat == 1),
  lo = as.numeric(cat == 2),
  hi = as.numeric(cat == 3)
)

test_that("normal outcome", {
  net <- set_ipd(dat, study = study, trt = trt, y = y)
  fit <- nma(net, regression = ~x1*.trt,
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10),
             prior_aux = half_normal(1))


  br2 <- bayes_R2(fit)

  expect_s3_class(br2, "nma_summary")
  expect_equal(br2$summary$parameter, "bayes_R2")

  # rstanarm::stan_glm(y ~ x1 * trt, data = dat) |> rstanarm::bayes_R2() |> summary()
  expect_equal(br2$summary$mean, c(bayes_R2 = 0.91), tolerance = 0.01)


  lr2 <- loo_R2(fit)

  expect_s3_class(lr2, "nma_summary")
  expect_equal(lr2$summary$parameter, "loo_R2")

  # rstanarm::stan_glm(y ~ x1 * trt, data = dat) |> rstanarm::loo_R2() |> summary()
  expect_equal(lr2$summary$mean, c(loo_R2 = 0.91), tolerance = 0.01)

  expect_equal(loo::loo_predictive_metric(fit),
               list(estimate = 0.2, se = 0.007),
               tolerance = 0.01)
})

test_that("poisson outcome", {
  net <- set_ipd(dat, study = study, trt = trt, r = rp, E = E)
  fit <- nma(net, regression = ~x1*.trt,
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10),
             prior_aux = half_normal(1))


  br2 <- bayes_R2(fit)

  expect_s3_class(br2, "nma_summary")
  expect_equal(br2$summary$parameter, "bayes_R2")

  # fitb <- brms::brm(r ~ x1 * trt + offset(log(E)), data = dat, family = "poisson")
  # brms::bayes_R2(fitb)
  expect_equal(br2$summary$mean, c(bayes_R2 = 0.74), tolerance = 0.01)


  lr2 <- loo_R2(fit)

  expect_s3_class(lr2, "nma_summary")
  expect_equal(lr2$summary$parameter, "loo_R2")

  # brms::loo_R2(fitb)
  expect_equal(lr2$summary$mean, c(loo_R2 = 0.70), tolerance = 0.01)
})


test_that("binary outcome", {
  net <- set_ipd(dat, study = study, trt = trt, r = r)
  fit <- nma(net, regression = ~x1*.trt,
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10))

  br2 <- bayes_R2(fit)

  expect_s3_class(br2, "nma_summary")
  expect_equal(br2$summary$parameter, "bayes_R2")

  # rstanarm::stan_glm(r ~ x1 * trt, data = dat, family = "binomial") |> rstanarm::bayes_R2() |> summary()
  expect_equal(br2$summary$mean, c(bayes_R2 = 0.12), tolerance = 0.01)


  lr2 <- suppressWarnings(loo_R2(fit))

  expect_s3_class(lr2, "nma_summary")
  expect_equal(lr2$summary$parameter, "loo_R2")

  # rstanarm::stan_glm(r ~ x1 * trt, data = dat, family = "binomial") |> rstanarm::loo_R2() |> summary()
  expect_equal(lr2$summary$mean, c(loo_R2 = 0.10), tolerance = 0.01)


  expect_equal(suppressWarnings(loo::loo_predictive_metric(fit)),
               list(estimate = 0.666, se = 0.021),
               tolerance = 0.01)
})

test_that("ordered outcome", {
  net <- set_ipd(dat, study = study, trt = trt, r = multi(r0, lo, hi))
  fit <- suppressWarnings(nma(net, regression = ~x1*.trt,
             prior_intercept = normal(0, 100),
             prior_trt = normal(0, 10),
             prior_reg = normal(0, 10),
             prior_aux = flat()))


  br2 <- bayes_R2(fit)

  expect_s3_class(br2, "nma_summary")
  expect_equal(br2$summary$parameter, "bayes_R2")

  # fitb <- brms::brm(cat ~ x1 * trt, data = dat, family = "cumulative")
  # brms::bayes_R2(fitb)  # note - assumes equidistance between categories
  #  0.68
  expect_equal(br2$summary$mean, c(bayes_R2 = 0.73), tolerance = 0.01)


  # lr2 <- loo_R2(fit)
  #
  # expect_s3_class(lr2, "nma_summary")
  # expect_equal(lr2$summary$parameter, "loo_R2")
  #
  # # brms::loo_R2(fitb)  # note - assumes equidistance between categories
  # #  0.67
  # expect_equal(lr2$summary$mean, c(loo_R2 = 0.67), tolerance = 0.01)
})
