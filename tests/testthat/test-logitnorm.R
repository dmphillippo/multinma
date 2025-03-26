skip_on_cran()

test_that("qlogitnorm", {
  p <- runif(30000)

  mu <- c(0.1, 1, 3)
  sigma <- c(0.8, 0.5, 1)
  q1 <- qlogitnorm(p, mu = mu, sigma = sigma)
  expect_true(all(q1 <= 1 & q1 >= 0))
  expect_equal(plogitnorm(q1, mu = mu, sigma = sigma), p)

  m <- c(0.2, 0.5, 0.8)
  s <- c(0.1, 0.2, 0.05)
  q2 <- matrix(qlogitnorm(p, mean = m, sd = s), nrow = 3)
  expect_true(all(q2 <= 1 & q2 >= 0))
  expect_equal(c(plogitnorm(q2, mean = m, sd = s)), p)

  expect_equal(rowMeans(q2), m, tol = 0.01)
  expect_equal(apply(q2, 1, sd), s, tol = 0.01)
})

test_that("plogitnorm", {
  q <- runif(1000, 0.05, 0.95)

  mu <- c(0.1, 1, 3)
  sigma <- c(0.8, 0.5, 1)
  p1 <- plogitnorm(q, mu = mu, sigma = sigma)
  expect_true(all(p1 <= 1 & p1 >= 0))
  expect_equal(qlogitnorm(p1, mu = mu, sigma = sigma), q)

  m <- c(0.2, 0.5, 0.8)
  s <- c(0.1, 0.2, 0.05)
  p2 <- plogitnorm(q, mean = m, sd = s)
  expect_true(all(p2 <= 1 & p2 >= 0))
  expect_equal(qlogitnorm(p2, mean = m, sd = s), q)
})
