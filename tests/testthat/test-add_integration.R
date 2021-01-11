library(multinma)
library(dplyr)
library(purrr)

test_that("expects nma_data object", {
  expect_error(add_integration("uh oh"), "No add_integration method defined")
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
         x3_mean = rnorm(1), x3_sd = runif(1, 0.5, 2)) %>%
  ungroup()

smkdummy_contrast <- smkdummy %>%  group_by(studyn) %>%
  mutate(arm = 1:n(), r = if_else(arm == 1, NA_real_, r)) %>%
  ungroup()

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

  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = "oh dear"), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 1.1), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = -5), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 0), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 1:2), m)

  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = "oh dear"), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 1.1), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = -5), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 0), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = 1:2), m)
})

test_that("int_args is named list", {
  m <- "should be a named list"

  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = "oh dear"), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = 1.1), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = NULL), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = list(a = 1, 2)), m)

  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = "oh dear"), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = 1.1), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = NULL), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), int_args = list(a = 1, 2)), m)
})

test_that("cor should be correlation matrix or NULL", {
  m <- "should be a correlation matrix or NULL"

  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = "a"), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = list()), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = 2), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = matrix(1:4)), m)
  expect_error(add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = matrix(1:4, nrow = 2)), m)

  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = "a"), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = list()), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = 2), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = matrix(1:4)), m)
  expect_error(add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbinom, 1, x2),
                               cor = matrix(1:4, nrow = 2)), m)
})

test_that("cor must be specified if no IPD", {
  expect_error(add_integration(smknet_agd,
                               x1 = distr(qnorm, mean = 1, sd = 1),
                               x2 = distr(qbinom, 1, x2),
                               cor = NULL),
               "Specify a correlation matrix")

  expect_error(add_integration(smkdummy,
                               x1 = distr(qnorm, mean = 1, sd = 1),
                               x2 = distr(qbinom, 1, x2),
                               cor = NULL),
               "Specify a correlation matrix")
})

test_that("cor is not required if only one covariate", {
  expect_s3_class(add_integration(smknet_agd, x1 = distr(qnorm, mean = 1, sd = 1), cor = NULL),
                  "mlnmr_data")

  s <- add_integration(smkdummy, x1 = distr(qnorm, mean = 1, sd = 1), cor = NULL)
  expect_s3_class(s, "data.frame")
  expect_true(rlang::has_name(s, ".int_x1"))
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

  expect_error(add_integration(smkdummy, cor = cormat, x1 = "a"), m)
  expect_error(add_integration(smkdummy, cor = cormat, x1 = list()), m)
  expect_error(add_integration(smkdummy, cor = cormat, "a"), m)
  expect_error(add_integration(smkdummy, cor = cormat,
                               x1 = distr(qnorm, mean = 1, sd = 1), x2 = list(), m))
  expect_error(add_integration(smkdummy, cor = cormat,
                               x1 = distr(qnorm, mean = 1, sd = 1), distr(qnorm, mean = 1, sd = 1), m))
  expect_error(add_integration(smkdummy, cor = cormat), paste("No covariate distributions specified.+", m))
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
  smknet_miss <- combine_network(smknet_ipd,
                                 set_agd_arm(smkdummy, studyn, trtn, r = r, n = n))
  smknet_miss$ipd[5, "x1"] <- NA
  expect_warning(add_integration(smknet_miss,
                                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                 x2 = distr(qbinom, size = 1, prob = x2)),
                 "Missing values.+complete cases")
})

test_that("replaces already present integration points with warning", {
  smknet_int <- add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd))
  m <- "Replacing integration points already present in network"

  expect_warning(add_integration(smknet_int, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd)), m)
  expect_warning(add_integration(smknet_int, x3 = distr(qnorm, mean = x3_mean, sd = x3_sd)), m)

  smkdummy_int <- add_integration(smkdummy, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd))
  m2 <- "Replacing integration points already present in data frame"

  expect_warning(add_integration(smkdummy_int, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd)), m2)
  expect_warning(add_integration(smkdummy_int, x3 = distr(qnorm, mean = x3_mean, sd = x3_sd)), m2)
})

make_sub <- function(x, n, sub) {
  x[sample.int(length(x), n)] <- sub
  return(x)
}

test_that("error if qfun produces NaN, NA, NULL, Inf", {
  m <- "Invalid integration points were generated"

  smkdummy_miss <- smkdummy %>% mutate(x1_mean = make_sub(x1_mean, 5, NA),
                                       x3_mean = make_sub(x3_mean, 5, NA))

  smkdummy_miss_contrast <- smkdummy_contrast %>%
    mutate(x1_mean = make_sub(x1_mean, 5, NA), x3_mean = make_sub(x3_mean, 5, NA))

  smknet_miss <- set_agd_arm(smkdummy_miss, studyn, trtn, r = r, n = n) %>%
    combine_network(smknet_ipd)
  smknet_missc <- set_agd_contrast(smkdummy_miss_contrast, studyn, trtn, y = r, se = n) %>%
    combine_network(smknet_ipd)

  expect_error(suppressWarnings(
    add_integration(smknet_miss,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qbinom, 1, x2),
                    x3 = distr(qnorm, x3_mean, x3_sd))),
               m)

  expect_error(suppressWarnings(
    add_integration(smknet_missc,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qbinom, 1, x2),
                    x3 = distr(qnorm, x3_mean, x3_sd))),
               m)

  smkdummy_bad <- smkdummy %>% mutate(x1_sd = make_sub(x1_sd, 5, -1),
                                      x3_sd = make_sub(x3_sd, 5, -1))

  smkdummy_bad_contrast <- smkdummy_contrast %>%
    mutate(x1_sd = make_sub(x1_sd, 5, -1), x3_sd = make_sub(x3_sd, 5, -1))

  smknet_bad <- set_agd_arm(smkdummy_bad, studyn, trtn, r = r, n = n) %>%
    combine_network(smknet_ipd)
  smknet_badc <- set_agd_contrast(smkdummy_bad_contrast, studyn, trtn, y = r, se = n) %>%
    combine_network(smknet_ipd)

  expect_error(suppressWarnings(
    add_integration(smknet_bad,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x3 = distr(qnorm, mean = x3_mean, sd = x3_sd))),
               m)

  expect_error(suppressWarnings(
    add_integration(smknet_badc,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x3 = distr(qnorm, mean = x3_mean, sd = x3_sd))),
               m)

})

test_that("integration point marginals and correlations are correct", {
  skip_on_cran()
  tol <- 0.005
  cor_tol <- 0.05
  n_int <- 10000

  # 1 covariate
  s1 <- add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = n_int)
  expect_equal(map_dbl(s1$agd_arm$.int_x1, mean), s1$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s1$agd_arm$.int_x1, sd), s1$agd_arm$x1_sd, tolerance = tol)

  # 2 covariates, IPD cor matrix
  s2 <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       n_int = n_int)
  expect_equal(map_dbl(s2$agd_arm$.int_x1, mean), s2$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s2$agd_arm$.int_x1, sd), s2$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s2$agd_arm$.int_x2, mean), s2$agd_arm$x2, tolerance = tol)
  # expect_equal(map2_dbl(s2$agd_arm$.int_x1, s2$agd_arm$.int_x3, cor, method = "spearman"),
  #              rep(x1_x2_cor, nrow(s2$agd_arm)), tolerance = cor_tol)
  expect_equal(map2_dbl(s2$agd_arm$.int_x1, s2$agd_arm$.int_x3, cor, method = "spearman"),
               rep(s2$int_cor[1, 3], nrow(s2$agd_arm)), tolerance = cor_tol)

  # 2 covariates, user cor matrix
  s3 <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       cor = cormat,
                       n_int = n_int)
  expect_equal(map_dbl(s3$agd_arm$.int_x1, mean), s3$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s3$agd_arm$.int_x1, sd), s3$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s3$agd_arm$.int_x2, mean), s3$agd_arm$x2, tolerance = tol)
  expect_equal(map2_dbl(s3$agd_arm$.int_x1, s3$agd_arm$.int_x3, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(s3$agd_arm)), tolerance = cor_tol)

  # Correlations between continuous and discrete seem hard to produce from the copula
  skip("Correlations between continuous and discrete covariates are difficult to recreate")
  # expect_equal(map2_dbl(s2$agd_arm$.int_x1, s2$agd_arm$.int_x2, cor, method = "spearman"),
  #              rep(x1_x2_cor, nrow(s2$agd_arm)), tolerance = cor_tol)
  expect_equal(map2_dbl(s2$agd_arm$.int_x1, s2$agd_arm$.int_x2, cor, method = "spearman"),
               rep(s2$int_cor[1, 2], nrow(s2$agd_arm)), tolerance = cor_tol)
  expect_equal(map2_dbl(s3$agd_arm$.int_x1, s3$agd_arm$.int_x2, cor, method = "spearman"),
               rep(x1_x2_cor, nrow(s3$agd_arm)), tolerance = cor_tol)
})
