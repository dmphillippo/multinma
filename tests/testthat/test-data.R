test_that("empty objects are produced", {
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

test_that("error if `data` does not inherit data.frame", {
  vec <- 1:5

  expect_error(set_ipd(vec), "should be a data frame")
  expect_error(set_agd_arm(vec), "should be a data frame")
  expect_error(set_agd_contrast(vec), "should be a data frame")
})
