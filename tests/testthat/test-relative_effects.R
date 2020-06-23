
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
  expect_error(relative_effects(smk_fit_RE, all_contrasts = "a"), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = 1), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = list()), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = NA), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = NULL), m)
})

test_that("trt_ref argument", {
  m <- "does not match a treatment in the network"
  expect_error(relative_effects(smk_fit_RE, trt_ref = "a"), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = 1), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = list("a")), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = NA), m)

  # expect_warning(relative_effects(smk_fit_RE, all_contrasts = TRUE, trt_ref = "Self-help"), "Ignoring `trt_ref`")
})

test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(relative_effects(smk_fit_RE, probs = "a"), m)
  expect_error(relative_effects(smk_fit_RE, probs = -1), m)
  expect_error(relative_effects(smk_fit_RE, probs = 1.5), m)
  expect_error(relative_effects(smk_fit_RE, probs = Inf), m)
  expect_error(relative_effects(smk_fit_RE, probs = list()), m)
  expect_error(relative_effects(smk_fit_RE, probs = NA), m)
  expect_error(relative_effects(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(relative_effects(smk_fit_RE, summary = "a"), m)
  expect_error(relative_effects(smk_fit_RE, summary = 1), m)
  expect_error(relative_effects(smk_fit_RE, summary = list()), m)
  expect_error(relative_effects(smk_fit_RE, summary = NA), m)
  expect_error(relative_effects(smk_fit_RE, summary = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(relative_effects(smk_fit_RE, newdata = "a"), m)
  expect_error(relative_effects(smk_fit_RE, newdata = 1), m)
  expect_error(relative_effects(smk_fit_RE, newdata = list()), m)
  expect_error(relative_effects(smk_fit_RE, newdata = NA), m)
})
