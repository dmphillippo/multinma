library(multinma)
library(dplyr)
library(purrr)

test_that("expects nma_data object", {
  expect_error(add_integration("uh oh"), "nma_data")
})

test_that("error on empty network", {
  expect_error(add_integration(set_ipd(data.frame())), "Empty network")
})

test_that("error if no AgD", {
  expect_error(add_integration(set_ipd(smoking, studyn, trtc, y = r)), "No aggregate data")
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
smknet <- combine_network(smknet_agd, smknet_ipd)

test_that("n_int should be a positive integer", {
  m <- "should be a positive integer"
  expect_error(add_integration(smknet, n_int = "oh dear"), m)
  expect_error(add_integration(smknet, n_int = 1.1), m)
  expect_error(add_integration(smknet, n_int = -5), m)
  expect_error(add_integration(smknet, n_int = 0), m)
  expect_error(add_integration(smknet, n_int = 1:2), m)
})

test_that("int_args is named list", {
  m <- "should be a named list"
  expect_error(add_integration(smknet, int_args = "oh dear"), m)
  expect_error(add_integration(smknet, int_args = 1.1), m)
  expect_error(add_integration(smknet, int_args = NULL), m)
  expect_error(add_integration(smknet, int_args = list(a = 1, 2)), m)
})

test_that("cor should be correlation matrix or NULL", {
  m <- "should be a correlation matrix or NULL"
  expect_error(add_integration(smknet, cor = "a"), m)
  expect_error(add_integration(smknet, cor = list()), m)
  expect_error(add_integration(smknet, cor = 2), m)
  expect_error(add_integration(smknet, cor = matrix(1:4)), m)
  expect_error(add_integration(smknet, cor = matrix(1:4, nrow = 2)), m)
})

test_that("cor must be specified if no IPD", {
  expect_error(add_integration(smknet_agd, x1 = distr(qnorm, mean = 1, sd = 1), cor = NULL),
               "Specify a correlation matrix")
})

test_that("covariate arguments should be named distr", {
  m <- "should be specified as named arguments using the function `distr`"
  expect_error(add_integration(smknet, cor = cormat, x1 = "a"), m)
  expect_error(add_integration(smknet, cor = cormat, x1 = list()), m)
  expect_error(add_integration(smknet, cor = cormat, "a"), m)
  expect_error(add_integration(smknet, cor = cormat,
                               x1 = distr(qnorm, mean = 1, sd = 1), x2 = list(), m))
  expect_error(add_integration(smknet, cor = cormat,
                               x1 = distr(qnorm, mean = 1, sd = 1), distr(qnorm, mean = 1, sd = 1), m))
  expect_error(add_integration(smknet, cor = cormat), paste("No covariate distributions specified.+", m))
})

test_that("covariate names must match IPD if provided", {
  m <- "Covariate name\\(s\\) not found in IPD: aaa"
  expect_error(add_integration(smknet, aaa = distr(qnorm, mean = x1_mean, sd = x1_sd)), m)
  expect_error(add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               aaa = distr(qbinom, size = 1, prob = x2)), m)
  expect_error(add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               aaa = distr(qbinom, size = 1, prob = x2),
                               cor = cormat), m)
})

test_that("warning if missing covariate values when calculating correlations from IPD", {
  smknet_miss <- smknet
  smknet_miss$ipd[5, "x1"] <- NA
  expect_warning(add_integration(smknet_miss,
                                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                 x2 = distr(qbinom, size = 1, prob = x2)),
                 "Missing values.+complete cases")
})

test_that("integration point marginals and correlations are correct", {
  skip_on_cran()
  tol <- 0.001
  cor_tol <- 0.05
  n_int <- 10000

  # 1 covariate
  a <- add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = n_int)
  expect_equal(map_dbl(a$agd_arm$.int_x1, mean), a$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(a$agd_arm$.int_x1, sd), a$agd_arm$x1_sd, tolerance = tol)

  # 2 covariates, IPD cor matrix
  b <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       n_int = n_int)
  expect_equal(map_dbl(b$agd_arm$.int_x1, mean), b$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(b$agd_arm$.int_x1, sd), b$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(b$agd_arm$.int_x2, mean), b$agd_arm$x2, tolerance = tol)
  expect_equal(map2_dbl(b$agd_arm$.int_x1, b$agd_arm$.int_x3, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(b$agd_arm)), tolerance = cor_tol)

  # 2 covariates, user cor matrix
  c <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       cor = cormat,
                       n_int = n_int)
  expect_equal(map_dbl(c$agd_arm$.int_x1, mean), c$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(c$agd_arm$.int_x1, sd), c$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(c$agd_arm$.int_x2, mean), c$agd_arm$x2, tolerance = tol)
  expect_equal(map2_dbl(c$agd_arm$.int_x1, c$agd_arm$.int_x3, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(c$agd_arm)), tolerance = cor_tol)

  # Correlations between continuous and discrete seem hard to produce from the copula
  skip("Correlations between continuous and discrete covariates are difficult to recreate")
  expect_equal(map2_dbl(b$agd_arm$.int_x1, b$agd_arm$.int_x2, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(b$agd_arm)), tolerance = cor_tol)
  expect_equal(map2_dbl(c$agd_arm$.int_x1, c$agd_arm$.int_x2, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(c$agd_arm)), tolerance = cor_tol)
})
