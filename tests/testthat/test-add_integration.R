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
         x2 = rbeta(1, 5, 10),
         x3_mean = rnorm(1), x3_sd = runif(1, 0.5, 2),
         x4 = rbeta(1, 10, 5)) %>%
  ungroup()

smkdummy_contrast <- smkdummy %>%
  group_by(studyn) %>%
  mutate(arm = 1:n(),
         cc = if_else(any(r == 0), 0.5, 0),
         r = r + cc,
         n = n + cc,
         lor = if_else(arm == 1, NA_real_, log(r * (first(n) - first(r)) / ((n - r) * first(r)))),
         se = sqrt(if_else(arm == 1, 1 / r + 1 / (n - r), 1 / r + 1 / (n - r) + 1 / first(r) + 1 / (first(n) - first(r))))) %>%
  ungroup()

ns_ipd <- 2

xall_cor <- 0.25
n_i <- 200
cormat <- matrix(xall_cor, nrow = 4, ncol = 4)
diag(cormat) <- 1

cop <- copula::normalCopula(copula::P2p(cormat), dim = 4, dispstr = "un")
u <- matrix(runif(n_i * 4 * ns_ipd), ncol = 4)
u_cor <- as.data.frame(copula::cCopula(u, cop, inverse = TRUE))

ipddummy <-
  tibble(studyn = rep(ns_agd, n_i * ns_ipd) + rep(1:ns_ipd, each = n_i),
         trtn = sample(1:2, n_i * ns_ipd, TRUE)) %>%
  mutate(x1 = qnorm(u_cor[,1]),
         x2 = qbinom(u_cor[,2], 1, 0.6),
         x3 = qnorm(u_cor[,3], 1, 0.5),
         x4 = qbinom(u_cor[,4], 1, 0.4),
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
  smknet_missc <- set_agd_contrast(smkdummy_miss_contrast, studyn, trtn, y = lor, se = se) %>%
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
  smknet_badc <- set_agd_contrast(smkdummy_bad_contrast, studyn, trtn, y = lor, se = se) %>%
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

test_that("correctly assess type of marginal distributions", {
  dtype <- get_distribution_type(
             x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
             x2 = distr(qgamma, mean = x3_mean, sd = x3_sd),
             x3 = distr(qbern, prob = x2),
             x4 = distr(qbinom, size = 1, prob = x2),
             x5 = distr(qbinom, size = n, prob = x2),
             x6 = distr(qbinom, 1, x2),
             x7 = distr(qpois, lambda = x3_sd),
             x8 = distr(function(p, ...) qnorm(p, ...), mean = x1_mean, sd = x1_sd),
             x9 = distr(function(p, ...) qpois(p, ...), lambda = 0.5),
             x10 = distr(function(p, ...) qbern(p, ...), prob = x2),
             data = smkdummy)
  expect_identical(dtype,
                   c(x1 = "continuous",
                     x2 = "continuous",
                     x3 = "binary",
                     x4 = "binary",
                     x5 = "discrete",
                     x6 = "binary",
                     x7 = "discrete",
                     x8 = "continuous",
                     x9 = "discrete",
                     x10 = "binary"))
})

test_that("cor_adjust logic is correct", {
  expect_error(add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "bad"),
               "`cor_adjust` must be one of")

  # Spearman should be the default for cor = NULL
  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2))$int_cor, "cor_adjust"),
                   "spearman")
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "spearman")$int_cor, "copula_cor"),
               cor_adjust_spearman(
                 add_integration(smknet,
                                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                 x2 = distr(qbern, x2),
                                 cor_adjust = "spearman")$int_cor,
                 types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "pearson")$int_cor, "cor_adjust"),
                   "pearson")
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "pearson")$int_cor, "copula_cor"),
               cor_adjust_pearson(
                 add_integration(smknet,
                                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                 x2 = distr(qbern, x2),
                                 cor_adjust = "pearson")$int_cor,
                 types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "legacy")$int_cor, "cor_adjust"),
                   "legacy")
  # Check that "legacy" does the right (wrong) thing (Spearman cor with no adjustment)
  expect_equal(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "legacy")$int_cor, "copula_cor"),
               add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "spearman")$int_cor,
               check.attributes = FALSE,
               tolerance = 0)
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "legacy")$int_cor, "copula_cor"),
               add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "legacy")$int_cor,
               check.attributes = FALSE,
               tolerance = 0)

  expect_error(add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "none"),
               'Cannot specify cor_adjust = "none"')

  # Pearson should be the default for user-provided cor
  cmat <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor = cmat)$int_cor, "cor_adjust"),
                   "pearson")
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor = cmat)$int_cor, "copula_cor"),
               cor_adjust_pearson(cmat, types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "spearman",
                                        cor = cmat)$int_cor, "cor_adjust"),
                   "spearman")
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "spearman",
                                    cor = cmat)$int_cor, "copula_cor"),
               cor_adjust_spearman(cmat, types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "none",
                                        cor = cmat)$int_cor, "cor_adjust"),
                   "none")
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "none",
                                    cor = cmat)$int_cor, "copula_cor"),
               add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "none",
                               cor = cmat)$int_cor,
               check.attributes = FALSE,
               tolerance = 0)
  expect_equal(add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "none",
                               cor = cmat)$int_cor,
               cmat,
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smknet,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "legacy",
                                        cor = cmat)$int_cor, "cor_adjust"),
                   "legacy")

  # Check that "legacy" does the right thing - same as "none" with user cor
  expect_equal(attr(add_integration(smknet,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "legacy",
                                    cor = cmat)$int_cor, "copula_cor"),
               add_integration(smknet,
                               x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                               x2 = distr(qbern, x2),
                               cor_adjust = "legacy",
                               cor = cmat)$int_cor,
               check.attributes = FALSE,
               tolerance = 0)

  # Pearson should be the default for user-provided cor
  expect_identical(attr(add_integration(smkdummy,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor = cmat), "cor_adjust"),
                   "pearson")
  expect_equal(attr(add_integration(smkdummy,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor = cmat), "copula_cor"),
               cor_adjust_pearson(cmat, types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smkdummy,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "spearman",
                                        cor = cmat), "cor_adjust"),
                   "spearman")
  expect_equal(attr(add_integration(smkdummy,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor = cmat,
                                    cor_adjust = "spearman"), "copula_cor"),
               cor_adjust_spearman(cmat, types = c("continuous", "binary")),
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smkdummy,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "none",
                                        cor = cmat), "cor_adjust"),
                   "none")
  expect_equal(attr(add_integration(smkdummy,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "none",
                                    cor = cmat), "copula_cor"),
               cmat,
               check.attributes = FALSE,
               tolerance = 0)

  expect_identical(attr(add_integration(smkdummy,
                                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                        x2 = distr(qbern, x2),
                                        cor_adjust = "legacy",
                                        cor = cmat), "cor_adjust"),
                   "legacy")

  # Check that "legacy" does the right thing - same as "none" with user cor
  expect_equal(attr(add_integration(smkdummy,
                                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                    x2 = distr(qbern, x2),
                                    cor_adjust = "legacy",
                                    cor = cmat), "copula_cor"),
               cmat,
               check.attributes = FALSE,
               tolerance = 0)

  # Check that reusing $int_cor from a network works automatically
  s1 <- add_integration(smknet,
                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                        x2 = distr(qbinom, size = 1, prob = x2),
                        x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                        x4 = distr(qbern, prob = x4))

  expect_identical(add_integration(smknet,
                                   x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                   x2 = distr(qbinom, size = 1, prob = x2),
                                   x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                                   x4 = distr(qbern, prob = x4),
                                   cor = s1$int_cor),
                   s1)

  s2 <- add_integration(smknet,
                        x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                        x2 = distr(qbinom, size = 1, prob = x2),
                        x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                        x4 = distr(qbern, prob = x4),
                        cor_adjust = "pearson")

  expect_identical(add_integration(smknet,
                                   x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                   x2 = distr(qbinom, size = 1, prob = x2),
                                   x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                                   x4 = distr(qbern, prob = x4),
                                   cor = s2$int_cor),
                   s2)
})

test_that("integration point marginals and correlations are correct", {
  skip_on_cran()
  tol <- 0.005
  cor_tol <- 0.05
  n_int <- 2^14

  # 1 covariate
  s1 <- add_integration(smknet, x1 = distr(qnorm, mean = x1_mean, sd = x1_sd), n_int = n_int)
  expect_equal(map_dbl(s1$agd_arm$.int_x1, mean), s1$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s1$agd_arm$.int_x1, sd), s1$agd_arm$x1_sd, tolerance = tol)

  # Multiple covariates, IPD cor matrix
  s2_spearman <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       x4 = distr(qbern, prob = x4),
                       n_int = n_int)

  # Mean, SD, proportion should match
  expect_equal(map_dbl(s2_spearman$agd_arm$.int_x1, mean), s2_spearman$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s2_spearman$agd_arm$.int_x1, sd), s2_spearman$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s2_spearman$agd_arm$.int_x2, mean), s2_spearman$agd_arm$x2, tolerance = tol)

  # Cont-cont correlations should match exactly for Spearman correlation
  expect_equal(map2_dbl(s2_spearman$agd_arm$.int_x1, s2_spearman$agd_arm$.int_x3, cor, method = "spearman"),
               rep(s2_spearman$int_cor[1, 3], nrow(s2_spearman$agd_arm)), tolerance = tol)

  # Other correlations should be approximate
  int_cors <- s2_spearman$agd_arm %>% rowwise() %>%
    mutate(int_cors = list(cor(cbind(.int_x1, .int_x2, .int_x3, .int_x4), method = "spearman"))) %>%
    pull(int_cors)
  # expect_equal(unlist(map(int_cors, ~.[upper.tri(.)])),
  #              rep(s2_spearman$int_cor[upper.tri(s2_spearman$int_cor)], times = nrow(s2_spearman$agd_arm)),
  #              tolerance = cor_tol)
  expect_lt(median(abs(unlist(map(int_cors, ~.[upper.tri(.)])) -
                   rep(s2_spearman$int_cor[upper.tri(s2_spearman$int_cor)], times = nrow(s2_spearman$agd_arm)))),
            cor_tol)

  s2_pearson <- add_integration(smknet,
                                 x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                 x2 = distr(qbinom, size = 1, prob = x2),
                                 x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                                x4 = distr(qbern, prob = x4),
                                cor_adjust = "pearson",
                                 n_int = n_int)
  expect_equal(map_dbl(s2_pearson$agd_arm$.int_x1, mean), s2_pearson$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s2_pearson$agd_arm$.int_x1, sd), s2_pearson$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s2_pearson$agd_arm$.int_x2, mean), s2_pearson$agd_arm$x2, tolerance = tol)

  # Normal-Normal correlations should match exactly for Pearson correlation
  expect_equal(map2_dbl(s2_pearson$agd_arm$.int_x1, s2_pearson$agd_arm$.int_x3, cor, method = "pearson"),
               rep(s2_pearson$int_cor[1, 3], nrow(s2_pearson$agd_arm)), tolerance = tol)

  # Other correlations should be approximate
  int_cors <- s2_pearson$agd_arm %>% rowwise() %>%
    mutate(int_cors = list(cor(cbind(.int_x1, .int_x2, .int_x3, .int_x4), method = "pearson"))) %>%
    pull(int_cors)
  # expect_equal(unlist(map(int_cors, ~.[upper.tri(.)])),
  #              rep(s2_pearson$int_cor[upper.tri(s2_pearson$int_cor)], times = nrow(s2_pearson$agd_arm)),
  #              tolerance = cor_tol)
  expect_lt(median(abs(unlist(map(int_cors, ~.[upper.tri(.)])) -
                         rep(s2_pearson$int_cor[upper.tri(s2_pearson$int_cor)], times = nrow(s2_pearson$agd_arm)))),
            cor_tol)

  # Multiple covariates, user cor matrix
  s3_pearson <- add_integration(smknet,
                       x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                       x2 = distr(qbinom, size = 1, prob = x2),
                       x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                       x4 = distr(qbern, prob = x4),
                       cor = cormat,
                       n_int = n_int)
  expect_equal(map_dbl(s3_pearson$agd_arm$.int_x1, mean), s3_pearson$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s3_pearson$agd_arm$.int_x1, sd), s3_pearson$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s3_pearson$agd_arm$.int_x2, mean), s3_pearson$agd_arm$x2, tolerance = tol)

  # Normal-Normal correlations should match exactly for Pearson correlation
  expect_equal(map2_dbl(s3_pearson$agd_arm$.int_x1, s3_pearson$agd_arm$.int_x3, cor, method = "pearson"),
               rep(cormat[1, 3], nrow(s3_pearson$agd_arm)), tolerance = tol)

  # Other correlations should be approximate
  int_cors <- s3_pearson$agd_arm %>% rowwise() %>%
    mutate(int_cors = list(cor(cbind(.int_x1, .int_x2, .int_x3, .int_x4), method = "pearson"))) %>%
    pull(int_cors)
  # expect_equal(unlist(map(int_cors, ~.[upper.tri(.)])),
  #              rep(s3_pearson$int_cor[upper.tri(s3_pearson$int_cor)], times = nrow(s3_pearson$agd_arm)),
  #              tolerance = cor_tol)
  expect_lt(median(abs(unlist(map(int_cors, ~.[upper.tri(.)])) -
                         rep(s3_pearson$int_cor[upper.tri(s3_pearson$int_cor)], times = nrow(s3_pearson$agd_arm)))),
            cor_tol)


  s3_spearman <- add_integration(smknet,
                                x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                                x2 = distr(qbinom, size = 1, prob = x2),
                                x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                                x4 = distr(qbern, prob = x4),
                                cor = cormat, cor_adjust = "spearman",
                                n_int = n_int)
  expect_equal(map_dbl(s3_spearman$agd_arm$.int_x1, mean), s3_spearman$agd_arm$x1_mean, tolerance = tol)
  expect_equal(map_dbl(s3_spearman$agd_arm$.int_x1, sd), s3_spearman$agd_arm$x1_sd, tolerance = tol)
  expect_equal(map_dbl(s3_spearman$agd_arm$.int_x2, mean), s3_spearman$agd_arm$x2, tolerance = tol)

  # Cont-cont correlations should match exactly for Spearman correlation
  expect_equal(map2_dbl(s3_spearman$agd_arm$.int_x1, s3_spearman$agd_arm$.int_x3, cor, method = "spearman"),
               rep(cormat[1, 3], nrow(s3_spearman$agd_arm)), tolerance = tol)

  # Other correlations should be approximate
  int_cors <- s3_spearman$agd_arm %>% rowwise() %>%
    mutate(int_cors = list(cor(cbind(.int_x1, .int_x2, .int_x3, .int_x4), method = "spearman"))) %>%
    pull(int_cors)
  # expect_equal(unlist(map(int_cors, ~.[upper.tri(.)])),
  #              rep(s3_spearman$int_cor[upper.tri(s3_spearman$int_cor)], times = nrow(s3_spearman$agd_arm)),
  #              tolerance = cor_tol)
  expect_lt(median(abs(unlist(map(int_cors, ~.[upper.tri(.)])) -
                         rep(s3_spearman$int_cor[upper.tri(s3_spearman$int_cor)], times = nrow(s3_spearman$agd_arm)))),
            cor_tol)
})

test_that("non positive definite cor matrices are fixed up", {
  badcor <- matrix(0.99, nrow = 4, ncol = 4)
  diag(badcor) <- 1

  expect_warning(
    add_integration(smknet,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qbinom, size = 1, prob = x2),
                    x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                    x4 = distr(qbern, prob = x4),
                    cor = badcor,
                    n_int = 10),
    "Adjusted correlation matrix not positive definite; using Matrix::nearPD()."
  )

  expect_warning(
    add_integration(smknet,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qbinom, size = 1, prob = x2),
                    x3 = distr(qnorm, mean = x3_mean, sd = x3_sd),
                    x4 = distr(qbern, prob = x4),
                    cor = badcor, cor_adjust = "spearman",
                    n_int = 10),
    "Adjusted correlation matrix not positive definite; using Matrix::nearPD()."
  )
})
