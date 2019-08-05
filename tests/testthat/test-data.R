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
