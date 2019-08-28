library(multinma)

test_that("expects nma_data object", {
  expect_error(add_integration("uh oh"), "nma_data")
})

test_that("error on empty network", {
  expect_error(add_integration(set_ipd(data.frame())), "Empty network")
})

smknet <- set_agd_arm(smoking, studyn, trtc, r = r, n = n)

test_that("n_int should be a positive integer", {
  m <- "should be a positive integer"
  expect_error(add_integration(smknet, n_int = "oh dear"), m)
  expect_error(add_integration(smknet, n_int = 1.1), m)
  expect_error(add_integration(smknet, n_int = -5), m)
  expect_error(add_integration(smknet, n_int = 0), m)
  expect_error(add_integration(smknet, n_int = 1:2), m)
})

test_that("int_args is named list", {
  m <- "should be a named list"
  expect_error(add_integration(smknet, int_args = "oh dear"), m)
  expect_error(add_integration(smknet, int_args = 1.1), m)
  expect_error(add_integration(smknet, int_args = NULL), m)
  expect_error(add_integration(smknet, int_args = list(a = 1, 2)), m)
})
})
