
test_that("check_prior_location gives suitable errors",{
  m <- "Prior location \\(mean\\) must be numeric"
  expect_error(check_prior_location("a"), m)
  expect_error(check_prior_location(list()), m)
  expect_error(check_prior_location(1:2), m)
  expect_error(check_prior_location("a", type = "median"), "Prior median must be numeric")
})

test_that("check_prior_scale gives suitable errors",{
  m <- "Prior scale \\(standard deviation\\) must be numeric"
  expect_error(check_prior_scale("a"), m)
  expect_error(check_prior_scale(list()), m)
  expect_error(check_prior_scale(1:2), m)
  expect_error(check_prior_scale("a", type = "variance"), "Prior variance must be numeric")
  expect_error(check_prior_scale(0), "Prior scale \\(standard deviation\\) must be strictly positive")
  expect_error(check_prior_scale(-1), "Prior scale \\(standard deviation\\) must be strictly positive")
})
