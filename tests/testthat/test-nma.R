library(multinma)
library(dplyr)

test_that("nma() wants a nma_data object", {
  m <- "Expecting an `nma_data` object"
  expect_error(nma(network = 1), m)
  expect_error(nma(network = list()), m)
  expect_error(nma(network = NULL), m)
})

smknet <- set_agd_arm(smoking, studyn, trtc, r = r, n = n)

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
  expect_error(nma(smknet_2, regression = y ~ x), "one-sided formula")
  expect_error(nma(smknet_2, regression = ~ a), "Failed to construct design matrix")
  expect_error(nma(smknet_2, regression = ~ a*.trt), "Failed to construct design matrix")
  expect_error(nma(smknet_2, regression = ~ (a + b)*.trt), "Failed to construct design matrix")
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
  m <- "missing values"
  expect_error(nma(smknet_agd_missx, regression = ~x1_mean), m)
  expect_error(nma(smknet_ipd_missx, regression = ~x1), m)
  expect_error(suppressWarnings(
    # Suppress warning about naive plug-in model
    nma(combine_network(smknet_agd_missx, smknet_ipd_missx),
        regression = ~x1)), m)
})
