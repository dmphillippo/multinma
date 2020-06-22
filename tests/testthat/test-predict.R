
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

test_that("baseline argument", {
  m <- "should be specified using distr"
  expect_error(predict(smk_fit_RE, baseline = "a"), m)
  expect_error(predict(smk_fit_RE, baseline = 1), m)
  expect_error(predict(smk_fit_RE, baseline = list("a")), m)
  expect_error(predict(smk_fit_RE, baseline = NA), m)
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
