library(multinma)

test_that("expects nma_data object", {
  expect_error(add_integration("uh oh"), "nma_data")
})
