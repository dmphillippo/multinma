test_that("combine_network error if not passed nma_data objects", {
  msg <- "Expecting to combine objects of class `nma_data`, created using set_* functions"

  expect_error(combine_network(1), msg)
  expect_error(combine_network(1, 2), msg)
  expect_error(combine_network(1, set_ipd(smoking[NA, ])), msg)
  expect_error(combine_network(set_ipd(smoking[NA, ]), 1), msg)
})
