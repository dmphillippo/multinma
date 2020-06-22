
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

p <- normal(0, 1)

test_that("summary.nma_prior digits argument", {
  m <- "must be a single integer"
  expect_error(summary(p, digits = "a"), m)
  expect_error(summary(p, digits = list()), m)
  expect_error(summary(p, digits = NA), m)
  expect_error(summary(p, digits = NULL), m)
  expect_error(summary(p, digits = 1:2), m)
  expect_error(summary(p, digits = 1.5), m)
  expect_error(summary(p, digits = Inf), m)
})

test_that("summary.nma_prior probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(summary(p, probs = "a"), m)
  expect_error(summary(p, probs = -1), m)
  expect_error(summary(p, probs = 1.5), m)
  expect_error(summary(p, probs = Inf), m)
  expect_error(summary(p, probs = list()), m)
  expect_error(summary(p, probs = NA), m)
  expect_error(summary(p, probs = NULL), m)
})

test_that("summary.nma_prior trunc argument", {
  m <- "length 2 numeric vector"
  expect_error(summary(p, trunc = "a"), m)
  expect_error(summary(p, trunc = -1), m)
  expect_error(summary(p, trunc = 1.5), m)
  expect_error(summary(p, trunc = Inf), m)
  expect_error(summary(p, trunc = list()), m)
  expect_error(summary(p, trunc = NA), m)
  expect_error(summary(p, trunc = 1:3), m)
  expect_error(summary(p, trunc = c(1, NA)), m)
})

test_that("summary.nma_prior sense check", {
  expect_output(summary(p, probs = 0.95, digits = 2), "95%.+between -1\\.96 and 1\\.96")
})
