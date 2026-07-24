library(dplyr)


test_that("baseline_synthesis() prior_intercept_sd must be valid", {
  m_baseline_sd <- "`prior_intercept_sd` must be"
  expect_error(baseline_synthesis(sa_net, prior_intercept_sd = 1), m_baseline_sd)
  expect_error(baseline_synthesis(sa_net, prior_intercept_sd = "a"), m_baseline_sd)
})

# test_that("baseline_synthesis() requires a disconnected network", {
#   expect_error(baseline_synthesis(sa_net), "`baseline_synthesis()` is only for disconnected networks.", fixed = TRUE)
# })

# Minimal disconnected, AgD-arm-only network: two studies, no shared treatment
disc_dat <- tibble(study = c("S1", "S1", "S2", "S2"),
                   trt   = c("A", "B", "C", "D"),
                   r     = c(50, 60, 55, 65),
                   n     = c(100, 100, 100, 100))
disc_net <- set_agd_arm(disc_dat, study, trt, r = r, n = n)

test_that("baseline_synthesis() summarises baseline parameters for the reference subnetwork only", {
  skip_on_cran()

  fit <- suppressWarnings(baseline_synthesis(disc_net, iter = 50))

  expect_s3_class(fit, "stan_baseline")
  expect_s3_class(fit, "stan_nma")

  params <- fit$baseline_summary$parameter
  expect_true(all(c("baseline_new", "baseline_mean", "baseline_sd") %in% params))
  # Only S1 (the reference treatment's subnetwork) contributes a mu[] row
  expect_equal(sum(grepl("^mu\\[", params)), 1)
  # d[] splits into one relative effect per subnetwork
  expect_equal(sum(grepl("^d\\[", params)), 2)
})

# Minimal disconnected, AgD-arm-only network with three subnetworks, no shared treatments
disc_dat_3 <- tibble(study = c("S1", "S1", "S2", "S2", "S3", "S3"),
                     trt   = c("A", "B", "C", "D", "E", "F"),
                     r     = c(50, 60, 55, 65, 48, 58),
                     n     = c(100, 100, 100, 100, 100, 100))
disc_net_3 <- set_agd_arm(disc_dat_3, study, trt, r = r, n = n)

test_that("baseline_synthesis() works with more than two subnetworks", {
  skip_on_cran()

  fit <- suppressWarnings(baseline_synthesis(disc_net_3, iter = 50))

  expect_s3_class(fit, "stan_baseline")

  params <- fit$baseline_summary$parameter
  # Only S1 (the reference treatment's subnetwork) contributes a mu[] row
  expect_equal(sum(grepl("^mu\\[", params)), 1)
  # d[] splits into one relative effect per subnetwork (S1 vs S2 vs S3)
  expect_equal(sum(grepl("^d\\[", params)), 3)
})
