library(multinma)

test_that("nma() wants a nma_data object", {
  m <- "Expecting an `nma_data` object"
  expect_error(nma(network = 1), m)
  expect_error(nma(network = list()), m)
  expect_error(nma(network = NULL), m)
})

smknet <- set_agd_arm(smoking, studyn, trtc, r = r, n = n)

test_that("nma() consistency argument must be one listed", {
  m <- "`consistency` must be"
  expect_error(nma(smknet, consistency = "con"), m)
  expect_error(nma(smknet, consistency = "a"), m)
  expect_error(nma(smknet, consistency = 1), m)
  expect_error(nma(smknet, consistency = NA), m)
  expect_error(nma(smknet, consistency = NULL), m)

  # Not needed - rlang::arg_match always returns a single string
  # expect_error(nma(smknet, consistency = c("consistency", "ume")), m)
})

test_that("nma() trt_effects argument must be one listed", {
  m <- "`trt_effects` must be"
  expect_error(nma(smknet, trt_effects = "fix"), m)
  expect_error(nma(smknet, trt_effects = "a"), m)
  expect_error(nma(smknet, trt_effects = 1), m)
  expect_error(nma(smknet, trt_effects = NA), m)
  expect_error(nma(smknet, trt_effects = NULL), m)

  # Not needed - rlang::arg_match always returns a single string
  # expect_error(nma(smknet, trt_effects = c("fixed", "random")), m)
})

