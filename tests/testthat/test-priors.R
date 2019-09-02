
test_that("check_prior_location gives suitable errors",{
  expect_error(check_prior_location("a"), "Prior mean must be numeric")
  expect_error(check_prior_location(list()), "Prior mean must be numeric")
  expect_error(check_prior_location(1:2), "Prior mean must be numeric, length 1")
  expect_error(check_prior_location("a", type = "median"), "Prior median must be numeric")
})

test_that("check_prior_scale gives suitable errors",{
  expect_error(check_prior_scale("a"), "Prior standard deviation must be numeric")
  expect_error(check_prior_scale(list()), "Prior standard deviation must be numeric")
  expect_error(check_prior_scale(1:2), "Prior standard deviation must be numeric, length 1")
  expect_error(check_prior_scale("a", type = "variance"), "Prior variance must be numeric")
  expect_error(check_prior_scale(0), "Prior standard deviation must be strictly positive")
  expect_error(check_prior_scale(-1), "Prior standard deviation must be strictly positive")
})
