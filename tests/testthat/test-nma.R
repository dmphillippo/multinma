library(multinma)
library(dplyr)

# Expose internal nma.fit()
nma.fit <- multinma:::nma.fit

test_that("nma() wants a nma_data object", {
  m <- "Expecting an `nma_data` object"
  expect_error(nma(network = 1), m)
  expect_error(nma(network = list()), m)
  expect_error(nma(network = NULL), m)
})

smknet <- set_agd_arm(smoking, studyn, trtc, r = r, n = n, trt_ref = "No intervention")

test_that("nma() consistency argument must be one listed", {
  m <- "`consistency` must be"
  expect_error(nma(smknet, consistency = "con"), m)
  expect_error(nma(smknet, consistency = "a"), m)
  expect_error(nma(smknet, consistency = 1), m)
  expect_error(nma(smknet, consistency = NA), m)
  expect_error(nma(smknet, consistency = NULL), m)

  # Not needed - rlang::arg_match always returns a single string
  # expect_error(nma(smknet, consistency = c("consistency", "ume")), m)
})

test_that("nma() trt_effects argument must be one listed", {
  m <- "`trt_effects` must be"
  expect_error(nma(smknet, trt_effects = "fix"), m)
  expect_error(nma(smknet, trt_effects = "a"), m)
  expect_error(nma(smknet, trt_effects = 1), m)
  expect_error(nma(smknet, trt_effects = NA), m)
  expect_error(nma(smknet, trt_effects = NULL), m)

  # Not needed - rlang::arg_match always returns a single string
  # expect_error(nma(smknet, trt_effects = c("fixed", "random")), m)
})

test_that("nma() likelihood must be valid", {
  m <- "`likelihood` should be"
  expect_error(nma(smknet, likelihood = 1), m)
  expect_error(nma(smknet, likelihood = "a"), m)
  expect_error(nma(smknet, likelihood = "norm"), m)
  expect_error(nma(smknet, likelihood = c("normal", "bernoulli")), m)
})

test_that("nma() link must be valid", {
  m <- "`link` should be"
  expect_error(nma(smknet, link = 1), m)
  expect_error(nma(smknet, link = "a"), m)
  expect_error(nma(smknet, link = "lo"), m)
  expect_error(nma(smknet, link = c("log", "identity")), m)
})

sa_net <- set_agd_contrast(social_anxiety, studyc, trtc, y = y, se = se, trt_class = classc, trt_ref = "Waitlist")

#Class effect parameter error checks
test_that("nma() class_effect must be valid", {
  m <- "`class_effects` must be"
  expect_error(nma(sa_net, class_effects = c("independent", "common")), m)
  expect_error(nma(sa_net, class_effects = 1), m)
  expect_error(nma(sa_net, class_effects = "random"), m)
})

test_that("nma() class_sd lists must be valid", {
  expect_error(nma(sa_net, class_effect = "exchangeable", class_sd = list(group1 = c("Exposure", "NonExistentClass"))), "Some classes listed in 'class_sd' are not present in the network", fixed = TRUE)
  expect_error(nma(sa_net, class_effect = "exchangeable", class_sd = list(group1 = c("Exposure", "NSSA"), group2 = c("NSSA", "SSRI/SNRI"))), "Some classes are listed in more than one shared standard deviation group in 'class_sd'")
})

test_that("nma() class_sd must be valid", {
  m_class_sd <- "`class_sd` must be"
  expect_error(nma(sa_net, class_effect = "exchangeable", class_sd = 1), m_class_sd)
  expect_error(nma(sa_net, class_effect = "exchangeable", class_sd = "a"), m_class_sd)
})

test_that("nma() gives an informative error if class_effect are specified without classes in the network", {
  m_no_classes <- "Setting `class_effects` requires treatment classes to be specified"
  expect_error(nma(smknet, class_effect = "common"), m_no_classes)
  expect_error(nma(smknet, class_effect = "exchangeable"), m_no_classes)
})

test_that("class_effects = 'common' updates network$treatments", {
  net_common <- nma(sa_net,
                    class_effects = "common",
                    trt_effects = "random",
                    prior_trt = normal(0, 100),
                    prior_het = half_normal(5),
                    test_grad = TRUE)

  # Check that net_common$network$treatments is updated to classes
  expect_equal(levels(net_common$network$treatments),
               levels(sa_net$classes),
               info = "network$treatments should be replaced by class factor levels")
})

# Make dummy covariate data for smoking network
ns_agd <- max(smoking$studyn)
smkdummy <-
  smoking %>%
  group_by(studyn) %>%
  mutate(x1_mean = rnorm(1), x1_sd = runif(1, 0.5, 2),
         x2 = runif(1),
         x3_mean = rnorm(1), x3_sd = runif(1, 0.5, 2))

ns_ipd <- 2

x1_x2_cor <- 0.4
n_i <- 200
cormat <- matrix(x1_x2_cor, nrow = 3, ncol = 3)
diag(cormat) <- 1

cop <- copula::normalCopula(copula::P2p(cormat), dim = 3, dispstr = "un")
u <- matrix(runif(n_i * 3 * ns_ipd), ncol = 3)
u_cor <- as.data.frame(copula::cCopula(u, cop, inverse = TRUE))

ipddummy <-
  tibble(studyn = rep(ns_agd, n_i * ns_ipd) + rep(1:ns_ipd, each = n_i),
         trtn = sample(1:2, n_i * ns_ipd, TRUE)) %>%
  mutate(x1 = qnorm(u_cor[,1]),
         x2 = qbinom(u_cor[,2], 1, 0.6),
         x3 = qnorm(u_cor[,3], 1, 0.5),
         r = rbinom(n(), 1, 0.2))

smknet_agd <- set_agd_arm(smkdummy, studyn, trtn, r = r, n = n)
smknet_ipd <- set_ipd(ipddummy, studyn, trtn, r = r)
smknet_2 <- combine_network(smknet_agd, smknet_ipd) %>%
  add_integration(x1 = distr(qnorm, x1_mean, x1_sd),
                  x2 = distr(qbinom, 1, x2),
                  x3 = distr(qnorm, x3_mean, x3_sd))

test_that("nma() regression formula is valid", {
  expect_error(nma(smknet_2, regression = y ~ x, center = FALSE), "one-sided formula")
  expect_error(nma(smknet_2, regression = ~ a, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ a*.trt, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ a^2*.trt, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ (a + b)*.trt, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variables "a" and "b" not found')

  expect_error(nma(smknet_2, regression = y ~ x, center = TRUE), "one-sided formula")
  expect_error(nma(smknet_2, regression = ~ a, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ a*.trt, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ a^2*.trt, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variable "a" not found')
  expect_error(nma(smknet_2, regression = ~ (a + b)*.trt, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), 'Regression variables "a" and "b" not found')
})


test_that("aux_regression ignored with warning", {
  surv_net <- set_ipd(ndmm_agd, study, trt, Surv = Surv(eventtime, status))

  expect_warning(f <- nma(smknet, aux_regression = ~.trt,
                     prior_intercept = normal(0, 10),
                     prior_trt = normal(0, 10),
                     prior_aux_reg = normal(0, 5), test_grad = TRUE),
                 "Ignoring `aux_regression`, only supported for survival")

  expect_null(f$aux_regression)

  expect_warning(f <- nma(surv_net, aux_regression = ~.trt,
                          likelihood = "exponential",
                          prior_intercept = normal(0, 10),
                          prior_trt = normal(0, 10), test_grad = TRUE),
                 "Ignoring `aux_regression`, no auxiliary parameters in exponential model")

  expect_null(f$aux_regression)
})

test_that("aux_by ignored with warning", {
  surv_net <- set_ipd(ndmm_agd, study, trt, Surv = Surv(eventtime, status))

  expect_warning(f <- nma(smknet, aux_by = ".trt",
                          prior_intercept = normal(0, 10),
                          prior_trt = normal(0, 10),
                          prior_aux_reg = normal(0, 5), test_grad = TRUE),
                 "Ignoring `aux_by`, only supported for survival")

  expect_null(f$aux_by)

  expect_warning(f <- nma(surv_net, aux_by = ".trt",
                          likelihood = "exponential",
                          prior_intercept = normal(0, 10),
                          prior_trt = normal(0, 10), test_grad = TRUE),
                 "Ignoring `aux_by`, no auxiliary parameters in exponential model")

  expect_null(f$aux_by)
})

make_na <- function(x, n) {
  x[sample.int(length(x), n)] <- NA
  return(x)
}

smknet_agd_missx <- smkdummy %>%
  mutate(x1_mean = make_na(x1_mean, 2), x1 = x1_mean) %>%
  set_agd_arm(studyn, trtn, r = r, n = n)

smknet_ipd_missx <- ipddummy %>%
  mutate(x1 = make_na(x1, 10)) %>%
  set_ipd(studyn, trtn, r = r)

test_that("nma() error if missing values in outcomes or predictors", {
  m <- "missing or infinite values"
  expect_error(nma(smknet_agd_missx, regression = ~x1_mean, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), m)
  expect_error(nma(smknet_agd_missx, regression = ~x1_mean, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), m)
  expect_error(nma(smknet_ipd_missx, regression = ~x1, center = FALSE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), m)
  expect_error(nma(smknet_ipd_missx, regression = ~x1, center = TRUE,
                   prior_intercept = normal(0, 10),
                   prior_trt = normal(0, 10),
                   prior_reg = normal(0, 5)), m)
  expect_error(suppressWarnings(
    # Suppress warning about naive plug-in model
    nma(combine_network(smknet_agd_missx, smknet_ipd_missx),
        regression = ~x1, center = FALSE)), m)
  expect_error(suppressWarnings(
    # Suppress warning about naive plug-in model
    nma(combine_network(smknet_agd_missx, smknet_ipd_missx),
        regression = ~x1, center = TRUE)), m)
})

smknet_3l <- combine_network(
  set_agd_arm(smkdummy %>% mutate(extra = 0.5, tclass = if_else(trtn == 1, 0, 1)),
              studyn, trtn, r = r, n = n, trt_class = tclass),
  set_ipd(ipddummy %>% mutate(extra = TRUE, tclass = if_else(trtn == 1, 0, 1)),
          studyn, trtn, r = r, trt_class = tclass),
  trt_ref = 1) %>%
  add_integration(x1 = distr(qnorm, x1_mean, x1_sd),
                  x2 = distr(qbinom, 1, x2),
                  x3 = distr(qnorm, x3_mean, x3_sd))

smknet_3c <- combine_network(
  set_agd_arm(smkdummy %>% mutate(extra = 0.5, tclass = if_else(trtn == 1, 0, 1)),
              studyn, trtn, r = r, n = n, trt_class = tclass),
  set_ipd(ipddummy %>% mutate(extra = "a", tclass = if_else(trtn == 1, 0, 1)),
          studyn, trtn, r = r, trt_class = tclass),
  trt_ref = 1) %>%
  add_integration(x1 = distr(qnorm, x1_mean, x1_sd),
                  x2 = distr(qbinom, 1, x2),
                  x3 = distr(qnorm, x3_mean, x3_sd))

test_that("nma() silently drops unnecessary columns", {
  expect_error(nma(smknet_3l, regression = ~(x1 + x2 + x3):.trt,
                    prior_intercept = normal(0, 10),
                    prior_trt = normal(0, 10),
                    prior_reg = normal(0, 5),
                    test_grad = TRUE), NA)

  expect_error(nma(smknet_3c, regression = ~(x1 + x2 + x3):.trt,
                    prior_intercept = normal(0, 10),
                    prior_trt = normal(0, 10),
                    prior_reg = normal(0, 5),
                    test_grad = TRUE), NA)
})

test_that("nma.fit() error if only one of x or y provided", {
  x <- matrix(1, nrow = 3, ncol = 3)
  colnames(x) <- c("a", "b", "c")

  y <- tibble(.y = 1:3)

  m <- "both be present or both NULL"
  expect_error(nma.fit(ipd_x = x), m)
  expect_error(nma.fit(ipd_y = y), m)
  expect_error(nma.fit(agd_arm_x = x), m)
  expect_error(nma.fit(agd_arm_y = y), m)
  expect_error(nma.fit(agd_contrast_x = x), "all be present or all NULL")
  expect_error(nma.fit(agd_contrast_y = y), "all be present or all NULL")
  expect_error(nma.fit(agd_contrast_Sigma = list()), "all be present or all NULL")
})

test_that("nma.fit() error if x and y dimension mismatch", {
  x <- matrix(1, nrow = 3, ncol = 3)
  colnames(x) <- c("a", "b", "c")

  y <- tibble(.y = 1:2)

  m <- "Number of rows.+do not match"
  expect_error(nma.fit(ipd_x = x, ipd_y = y), m)
  expect_error(nma.fit(agd_arm_x = x, agd_arm_y = y, n_int = 1), m)
  expect_error(nma.fit(agd_contrast_x = x, agd_contrast_y = y, agd_contrast_Sigma = list(), n_int = 1), m)
})

test_that("nma.fit() error if x column names different", {
  x1 <- x2 <- matrix(1, nrow = 3, ncol = 3)
  colnames(x1) <- c(".study1", ".trt1", "c")
  colnames(x2) <- c(".study1", ".trt1", "D")
  x3 <- x1[, 1:2]

  y <- tibble(.y = 1:3, .se = 1:3)

  Sigma <- list(matrix(1, 3, 3))

  m <- "Non-matching columns"
  expect_error(nma.fit(ipd_x = x1, ipd_y = y,
                       agd_arm_x = x2, agd_arm_y = y,
                       n_int = 1), m)
  expect_error(nma.fit(ipd_x = x1, ipd_y = y,
                       agd_contrast_x = x2, agd_contrast_y = y, agd_contrast_Sigma = Sigma,
                       n_int = 1), m)
  expect_error(nma.fit(agd_arm_x = x1, agd_arm_y = y,
                       agd_contrast_x = x2, agd_contrast_y = y, agd_contrast_Sigma = Sigma,
                       n_int = 1), m)
  expect_error(nma.fit(ipd_x = x1, ipd_y = y,
                       agd_arm_x = x3, agd_arm_y = y,
                       n_int = 1), m)
  expect_error(nma.fit(ipd_x = x1, ipd_y = y,
                       agd_contrast_x = x3, agd_contrast_y = y, agd_contrast_Sigma = Sigma,
                       n_int = 1), m)
  expect_error(nma.fit(agd_arm_x = x2, agd_arm_y = y,
                       agd_contrast_x = x3, agd_contrast_y = y, agd_contrast_Sigma = Sigma,
                       n_int = 1), m)
})

test_that("nma.fit() error if agd_contrast_Sigma is not right dimensions", {
  x1 <- matrix(1, nrow = 3, ncol = 3)
  colnames(x1) <- c(".study1", ".trt1", "x")

  y <- tibble(.y = 1:3, .se = 1:3)

  Sigma1 <- list(matrix(1, 3, 3), matrix(1, 1, 1))
  Sigma2 <- list(matrix(1, 4, 4))

  expect_error(nma.fit(agd_contrast_x = x1, agd_contrast_y = y,
                       agd_contrast_Sigma = Sigma1, likelihood = "normal", link = "identity",
                       prior_intercept = normal(0, 10),
                       prior_trt = normal(0, 10),
                       prior_reg = normal(0, 5),
                       prior_het = normal(0, 1),
                       n_int = 1), "Dimensions of `agd_contrast_Sigma`.+do not match")
  expect_error(nma.fit(agd_contrast_x = x1, agd_contrast_y = y,
                       agd_contrast_Sigma = Sigma2, likelihood = "normal", link = "identity",
                       prior_intercept = normal(0, 10),
                       prior_trt = normal(0, 10),
                       prior_reg = normal(0, 5),
                       prior_het = normal(0, 1),
                       n_int = 1), "Dimensions of `agd_contrast_Sigma`.+do not match")
})


single_study_a <- tibble(study = "A", trt = c("a", "b"), r = 500, n = 1000)
net_ss_a <- set_agd_arm(single_study_a, study, trt, r = r, n = n)

single_study_c <- tibble(study = "A", trt = c("a", "b"), y = c(NA, 0), se = c(NA, sqrt(4/500)))
net_ss_c <- set_agd_contrast(single_study_c, study, trt, y = y, se = se)

single_study_i <- tibble(study = "A", trt = rep(c("a", "b"), each = 100), y = c(rnorm(100, 0, 0.1), rnorm(100, 1, 0.1)))
net_ss_i <- set_ipd(single_study_i, study, trt, y = y)

test_that("nma() doesn't fail with a single study", {
  skip_on_cran()

  fit_ss_a <- nma(net_ss_a, prior_intercept = normal(0, 10), prior_trt = normal(0, 10))
  expect_s3_class(fit_ss_a, "stan_nma")
  d_ss_a <- relative_effects(fit_ss_a)
  expect_equivalent(d_ss_a$summary$mean, 0, tol = 0.01)
  expect_equivalent(d_ss_a$summary$sd, sqrt(4/500), tol = 0.01)

  fit_ss_c <- nma(net_ss_c, prior_intercept = normal(0, 10), prior_trt = normal(0, 10))
  expect_s3_class(fit_ss_c, "stan_nma")
  d_ss_c <- relative_effects(fit_ss_c)
  expect_equivalent(d_ss_c$summary$mean, 0, tol = 0.01)
  expect_equivalent(d_ss_c$summary$sd, sqrt(4/500), tol = 0.01)

  fit_ss_i <- nma(net_ss_i, prior_intercept = normal(0, 10), prior_trt = normal(0, 10), prior_aux = half_normal(2))
  expect_s3_class(fit_ss_i, "stan_nma")
  d_ss_i <- relative_effects(fit_ss_i)
  expect_equivalent(d_ss_i$summary$mean, 1, tol = 0.1)
  expect_equivalent(d_ss_i$summary$sd, sqrt(0.1^2 + 0.1^2)/sqrt(200), tol = 0.01)
})

test_that("nma() gives warnings for default priors", {
  m <- "Prior distributions were left at default values:"

  # Use test_grad = TRUE to avoid running any samples, just check logic works

  expect_warning(nma(smknet, test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+"))
  expect_warning(nma(smknet, prior_intercept = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_trt.+"))
  expect_warning(nma(smknet, prior_trt = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_intercept.+"))

  expect_warning(nma(smknet, trt_effects = "random", test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+", "prior_het.+"))
  expect_warning(nma(smknet, trt_effects = "random", prior_intercept = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_trt.+", "prior_het.+"))
  expect_warning(nma(smknet, trt_effects = "random", prior_trt = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_het.+"))
  expect_warning(nma(smknet, trt_effects = "random", prior_het = half_normal(1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+"))

  smknet_yi <- set_ipd(smoking, studyn, trtc, y = r)

  expect_warning(nma(smknet_yi, test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, prior_intercept = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_trt.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, prior_trt = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, prior_aux = half_normal(1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+"))

  expect_warning(nma(smknet_yi, trt_effects = "random", test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+", "prior_het.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, trt_effects = "random", prior_intercept = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_trt.+", "prior_het.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, trt_effects = "random", prior_trt = normal(0, 1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_het.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, trt_effects = "random", prior_het = half_normal(1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+", "prior_aux.+"))
  expect_warning(nma(smknet_yi, trt_effects = "random", prior_aux = half_normal(1), test_grad = TRUE), paste0(m, ".+prior_intercept.+", "prior_trt.+", "prior_het.+"))

})

test_that("nma() error with incompatible priors", {
  m <- "Invalid `prior_.+`\\. Suitable distributions are"

  expect_error(nma(smknet, prior_intercept = half_normal(1), prior_trt = normal(0, 1)), m)
  expect_error(nma(smknet, prior_intercept = normal(0, 1), prior_trt = half_normal(1)), m)
})

test_that("nma() prior checks work", {
  expect_error(nma(smknet, prior_intercept = "bad"),
               "`prior_intercept` must be a prior distribution")
  expect_error(nma(smknet, prior_intercept = list(normal(0, 1))),
               "`prior_intercept` must be a prior distribution")

  ndmm_net <- set_ipd(ndmm_ipd, studyf, trtf, Surv = Surv(eventtime, status))

  expect_error(nma(ndmm_net,
                   likelihood = "gengamma",
                   prior_intercept = normal(0, 100),
                   prior_trt = normal(0, 100),
                   prior_aux = letters),
               "`prior_aux` must be a named list of prior distributions")

  expect_error(nma(ndmm_net,
                   likelihood = "gengamma",
                   prior_intercept = normal(0, 100),
                   prior_trt = normal(0, 100),
                   prior_aux = list(a = normal(0, 1))),
               "Missing named elements with priors for sigma and k")

  expect_error(nma(ndmm_net,
                   likelihood = "gengamma",
                   prior_intercept = normal(0, 100),
                   prior_trt = normal(0, 100),
                   prior_aux = list(k = normal(0, 1))),
               "Missing named elements with priors for sigma")
})

test_that("rstan R-hat and ESS warnings are captured correctly", {
  # Need to capture these for numerical integration checks, but fragile since warning text may change
  local_edition(3)

  expect_warning(
  expect_warning(
  expect_warning(
    nma(smknet,
        prior_intercept = normal(0, 100),
        prior_trt = normal(0, 100),
        iter = 100,
        seed = 1234,
        init = 0),
    "The largest R-hat is"),
    "Bulk Effective Sample(s?) Size \\(ESS\\) is too low"),
    "Tail Effective Sample(s?) Size \\(ESS\\) is too low")
})
