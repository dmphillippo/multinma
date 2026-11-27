test_that("[dpq]gamma() respect the scale argument", {
  x <- c(0.5, 2, 5)
  p <- c(0.1, 0.5, 0.9)

  expect_equal(qgamma(p, shape = 2, scale = 3), stats::qgamma(p, shape = 2, scale = 3))
  expect_equal(pgamma(x, shape = 2, scale = 3), stats::pgamma(x, shape = 2, scale = 3))
  expect_equal(dgamma(x, shape = 2, scale = 3), stats::dgamma(x, shape = 2, scale = 3))

  expect_equal(qgamma(p, shape = 2, rate = 0.5), stats::qgamma(p, shape = 2, rate = 0.5))
  expect_equal(pgamma(x, shape = 2, rate = 0.5), stats::pgamma(x, shape = 2, rate = 0.5))
  expect_equal(dgamma(x, shape = 2, rate = 0.5), stats::dgamma(x, shape = 2, rate = 0.5))
})
