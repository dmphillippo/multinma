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

test_that("cor should be correlation matrix or NULL", {
  m <- "should be a correlation matrix or NULL"
  expect_error(add_integration(smknet, cor = "a"), m)
  expect_error(add_integration(smknet, cor = list()), m)
  expect_error(add_integration(smknet, cor = 2), m)
  expect_error(add_integration(smknet, cor = matrix(1:4)), m)
  expect_error(add_integration(smknet, cor = matrix(1:4, nrow = 2)), m)
})

test_that("covariate arguments should be named distr", {
  m <- "should be specified as named arguments using the function `distr`"
  expect_error(add_integration(smknet, x1 = "a"), m)
  expect_error(add_integration(smknet, x1 = list()), m)
  expect_error(add_integration(smknet, "a"), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = 1, sd = 1), x2 = list(), m))
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = 1, sd = 1), distr(qnorm, mean = 1, sd = 1), m))
})
