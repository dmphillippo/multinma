test_that("set_* produces empty nma_data objects", {
  empty_nma_data <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = NULL,
         studies = NULL), class = "nma_data")

  expect_equal(set_ipd(smoking[NA, ]), empty_nma_data)
  expect_equal(set_agd_arm(smoking[NA, ]), empty_nma_data)
  expect_equal(set_agd_contrast(smoking[NA, ]), empty_nma_data)
})

test_that("set_* error if data does not inherit data.frame", {
  vec <- 1:5
  msg <- "Argument `data` should be a data frame"

  expect_error(set_ipd(vec), msg)
  expect_error(set_agd_arm(vec), msg)
  expect_error(set_agd_contrast(vec), msg)
})

test_that("combine_network error if not passed nma_data objects", {
  msg <- "Expecting to combine objects of class `nma_data`, created using set_* functions"

  expect_error(combine_network(1), msg)
  expect_error(combine_network(1, 2), msg)
  expect_error(combine_network(1, set_ipd(smoking[NA, ])), msg)
  expect_error(combine_network(set_ipd(smoking[NA, ]), 1), msg)
})
