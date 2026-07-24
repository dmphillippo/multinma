library(dplyr)


test_that("baseline_synthesis() prior_intercept_sd must be valid", {
  m_baseline_sd <- "`prior_intercept_sd` must be"
  expect_error(baseline_synthesis(sa_net, prior_intercept_sd = 1), m_baseline_sd)
  expect_error(baseline_synthesis(sa_net, prior_intercept_sd = "a"), m_baseline_sd)
})

# Minimal disconnected, AgD-arm-only network: two studies, no shared treatment
disc_dat <- tibble(study = c("S1", "S1", "S2", "S2"),
                   trt   = c("A", "B", "C", "D"),
                   r     = c(50, 60, 55, 65),
                   n     = c(100, 100, 100, 100))
disc_net <- set_agd_arm(disc_dat, study, trt, r = r, n = n)

test_that("correct class", {
  fit <- suppressWarnings(baseline_synthesis(disc_net, test_grad = TRUE))
  expect_s3_class(fit, "stan_baseline")
  expect_s3_class(fit, "stan_nma")
})


skip_on_cran()


test_that("default summmary for correct parameters", {
  fit <- suppressWarnings(baseline_synthesis(disc_net, iter = 10))
  s <- as.data.frame(summary(fit))

  expect_true(all(c("baseline_new", "baseline_mean", "baseline_sd") %in% s$parameter))

  skip("To re-implement")

  # Only S1 (the reference treatment's subnetwork) contributes a mu[] row
  expect_equal(sum(grepl("^mu\\[", s$parameter)), 1)
  # d[] splits into one relative effect per subnetwork
  expect_equal(sum(grepl("^d\\[", s$parameter)), 2)
})

# Minimal disconnected, AgD-arm-only network with three subnetworks, no shared treatments
disc_dat_3 <- tibble(study = c("S1", "S1", "S2", "S2", "S3", "S3"),
                     trt   = c("A", "B", "C", "D", "E", "F"),
                     r     = c(50, 60, 55, 65, 48, 58),
                     n     = c(100, 100, 100, 100, 100, 100))
disc_net_3 <- set_agd_arm(disc_dat_3, study, trt, r = r, n = n)

test_that("baseline_synthesis() works with more than two subnetworks", {
  fit <- suppressWarnings(baseline_synthesis(disc_net_3, iter = 10))

  s <- as.data.frame(summary(fit))

  expect_true(all(c("baseline_new", "baseline_mean", "baseline_sd") %in% s$parameter))

  skip("To re-implement")

  # Only S1 (the reference treatment's subnetwork) contributes a mu[] row
  expect_equal(sum(grepl("^mu\\[", s$parameter)), 1)
  # d[] splits into one relative effect per subnetwork (S1 vs S2 vs S3)
  expect_equal(sum(grepl("^d\\[", s$parameter)), 3)
})

# Check IPD + AgD agree
bdat_a <- data.frame(study = c("S1", "S1", "S2", "S3", "S4", "S5"),
                     trt   = c("A", "B", "A", "A", "A", "C"),
                     r     = c(5, 6, 1, 9, 4, 7),
                     n     = c(10, 10, 10, 10, 10, 10))

bnet_a <- set_agd_arm(bdat_a, study, trt, r = r, n = n,
                      allow_singlearm_studies = TRUE)

bdat_i <- rowwise(bdat_a) %>%
  mutate(ri = list(c(rep(1, r), rep(0, n - r)))) %>%
  tidyr::unnest(cols = "ri")

bnet_i <- set_ipd(bdat_i, study, trt, r = ri,
                  allow_singlearm_studies = TRUE)

bnet_ai <- combine_network(
  set_agd_arm(bdat_a %>% filter(study == "S1"), study, trt, r = r, n = n,
              allow_singlearm_studies = TRUE),
  set_ipd(bdat_i %>% filter(study != "S1"), study, trt, r = ri,
          allow_singlearm_studies = TRUE)
)

test_that("AgD, IPD, and mixed analysis identical", {

  fit_a <- suppressWarnings(
    baseline_synthesis(bnet_a,
                       prior_intercept = normal(0, 10),
                       prior_intercept_sd = half_normal(0.1),
                       prior_trt = normal(0, 1),
                       iter = 10000))

  fit_i <- suppressWarnings(
             baseline_synthesis(bnet_i,
                              prior_intercept = normal(0, 10),
                              prior_intercept_sd = half_normal(0.1),
                              prior_trt = normal(0, 1),
                              iter = 10000))

  fit_ai <- suppressWarnings(
             baseline_synthesis(bnet_ai,
                              prior_intercept = normal(0, 10),
                              prior_intercept_sd = half_normal(0.1),
                              prior_trt = normal(0, 1),
                              iter = 10000))

  s_a <- as.data.frame(summary(fit_a)) %>% select(-"Bulk_ESS", -"Tail_ESS", -"Rhat")
  s_i <- as.data.frame(summary(fit_i)) %>% select(-"Bulk_ESS", -"Tail_ESS", -"Rhat")
  s_ai <- as.data.frame(summary(fit_ai)) %>% select(-"Bulk_ESS", -"Tail_ESS", -"Rhat")

  expect_equal(s_a, s_i, tolerance = 0.05)
  expect_equal(s_a, s_ai, tolerance = 0.05)

})

# Check synthesis with Normal likelihood
test_that("correct posterior - normal likelihood, disconnected", {

  cdat <- data.frame(study = c("S1", "S1", "S2", "S2", "S3", "S4", "S5", "S5"),
                     trt   = c("A", "B", "A", "D", "A", "E", "F", "G"),
                     y = c(1, 2, 1.1, 3, 0.8, 0.5, 0.6, 0.7), # runif(8, 0.5, 1.5),
                     se = runif(8, 0.1, 0.25))

  cnet <- set_agd_arm(cdat, study, trt, y = y, se = se,
                      allow_singlearm_studies = TRUE)

  cfit <- suppressWarnings(baseline_synthesis(cnet,
                             prior_intercept = normal(0, 10),
                             prior_intercept_sd = half_normal(0.1),
                             prior_trt = normal(0, 1),
                             iter = 10000))

  bl <- filter(cdat, trt == "A")
  tau_mu <- as.data.frame(summary(cfit, pars = "baseline_sd"))$mean

  # baseline mean
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_mean"))$mean,
               weighted.mean(bl$y, 1 / (bl$se^2 + tau_mu^2)),
               tolerance = 0.05)
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_mean"))$sd,
               sqrt(1/sum(1 / (bl$se^2 + tau_mu^2))),
               tolerance = 0.05)

  # predictive dist
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_new"))$mean,
               as.data.frame(summary(cfit, pars = "baseline_mean"))$mean,
               tolerance = 0.05)
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_new"))$sd,
               as.data.frame(summary(cfit, pars = "baseline_mean"))$sd +
                 as.data.frame(summary(cfit, pars = "baseline_sd"))$mean,
               tolerance = 0.05)

})

test_that("correct posterior - normal likelihood, connected", {

  cdat <- data.frame(study = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
                     trt   = c("A", "B", "A", "D", "A", "E", "F", "G"),
                     y = c(1, 2, 1.1, 3, 0.8, 0.5, 0.6, 0.7), # runif(8, 0.5, 1.5),
                     se = runif(8, 0.1, 0.25)) %>% filter(study %in% c("S1", "S2", "S3"))

  cnet <- set_agd_arm(cdat, study, trt, y = y, se = se)

  cfit <- suppressWarnings(baseline_synthesis(cnet,
                             prior_intercept = normal(0, 10),
                             prior_intercept_sd = half_normal(0.1),
                             prior_trt = normal(0, 1),
                             iter = 10000))

  bl <- filter(cdat, trt == "A")
  tau_mu <- as.data.frame(summary(cfit, pars = "baseline_sd"))$mean

  # baseline mean
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_mean"))$mean,
               weighted.mean(bl$y, 1 / (bl$se^2 + tau_mu^2)),
               tolerance = 0.05)
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_mean"))$sd,
               sqrt(1/sum(1 / (bl$se^2 + tau_mu^2))),
               tolerance = 0.05)

  # predictive dist
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_new"))$mean,
               as.data.frame(summary(cfit, pars = "baseline_mean"))$mean,
               tolerance = 0.05)
  expect_equal(as.data.frame(summary(cfit, pars = "baseline_new"))$sd,
               as.data.frame(summary(cfit, pars = "baseline_mean"))$sd +
                 as.data.frame(summary(cfit, pars = "baseline_sd"))$mean,
               tolerance = 0.05)

})

test_that("TSD5 smoking cessation - simultaneous modelling", {
  smknet <- set_agd_arm(smoking, studyn, trtc, r = r, n = n,
                        trt_ref = "No intervention")
  fit <- baseline_synthesis(smknet,
                            trt_effects = "random",
                            prior_intercept = normal(scale = 100),
                            prior_trt = normal(scale = 100),
                            prior_het = normal(scale = 5),
                            prior_intercept_sd = half_normal(2.5))


  tol <- 0.05

  s_mean <- as.data.frame(summary(fit, pars = "baseline_mean"))
  expect_equal(s_mean$mean, -2.49, tolerance = tol)
  expect_equal(s_mean$sd, 0.13, tolerance = tol)
  expect_equal(s_mean$`2.5%`, -2.75, tolerance = tol)
  expect_equal(s_mean$`97.5%`, -2.25, tolerance = tol)

  s_sd <- as.data.frame(summary(fit, pars = "baseline_sd"))
  expect_equal(s_sd$`50%`, 0.45, tolerance = tol)
  expect_equal(s_sd$sd, 0.11, tolerance = tol)
  expect_equal(s_sd$`2.5%`, 0.29, tolerance = tol)
  expect_equal(s_sd$`97.5%`, 0.71, tolerance = tol)

  s_new <- as.data.frame(summary(fit, pars = "baseline_new"))
  expect_equal(s_new$mean, -2.49, tolerance = tol)
  expect_equal(s_new$sd, 0.49, tolerance = tol)
  expect_equal(s_new$`2.5%`, -3.48, tolerance = tol)
  expect_equal(s_new$`97.5%`, -1.52, tolerance = tol)
})


test_that("TSD5 smoking cessation - separate modelling", {
  smknet <- set_agd_arm(smoking %>% filter(trtc == "No intervention"),
                        studyn, trtc, r = r, n = n,
                        allow_singlearm_studies = TRUE)
  fit <- baseline_synthesis(smknet,
                            trt_effects = "random",
                            prior_intercept = normal(scale = 100),
                            prior_trt = normal(scale = 100),
                            prior_het = normal(scale = 5),
                            prior_intercept_sd = half_normal(2.5))


  tol <- 0.05

  s_mean <- as.data.frame(summary(fit, pars = "baseline_mean"))
  expect_equal(s_mean$mean, -2.59, tolerance = tol)
  expect_equal(s_mean$sd, 0.16, tolerance = tol)
  expect_equal(s_mean$`2.5%`, -2.94, tolerance = tol)
  expect_equal(s_mean$`97.5%`, -2.30, tolerance = tol)

  s_sd <- as.data.frame(summary(fit, pars = "baseline_sd"))
  expect_equal(s_sd$`50%`, 0.54, tolerance = tol)
  expect_equal(s_sd$sd, 0.16, tolerance = tol)
  expect_equal(s_sd$`2.5%`, 0.32, tolerance = tol)
  expect_equal(s_sd$`97.5%`, 0.93, tolerance = tol)

  s_new <- as.data.frame(summary(fit, pars = "baseline_new"))
  expect_equal(s_new$mean, -2.59, tolerance = tol)
  expect_equal(s_new$sd, 0.60, tolerance = tol)
  expect_equal(s_new$`2.5%`, -3.82, tolerance = tol)
  expect_equal(s_new$`97.5%`, -1.41, tolerance = tol)
})
