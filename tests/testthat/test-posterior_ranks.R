
smk_net <- set_agd_arm(smoking,
                       study = studyn,
                       trt = trtc,
                       r = r,
                       n = n,
                       trt_ref = "No intervention")

# Only test gradients, no sampling
smk_fit_RE <- nma(smk_net,
                  trt_effects = "random",
                  prior_intercept = normal(scale = 100),
                  prior_trt = normal(scale = 100),
                  prior_het = normal(scale = 5),
                  test_grad = TRUE)


test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(posterior_ranks(smk_fit_RE, probs = "a"), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = -1), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = 1.5), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = Inf), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = list()), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = NA), m)
  expect_error(posterior_ranks(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_ranks(smk_fit_RE, summary = "a"), m)
  expect_error(posterior_ranks(smk_fit_RE, summary = 1), m)
  expect_error(posterior_ranks(smk_fit_RE, summary = list()), m)
  expect_error(posterior_ranks(smk_fit_RE, summary = NA), m)
  expect_error(posterior_ranks(smk_fit_RE, summary = NULL), m)
})

test_that("sucra argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_ranks(smk_fit_RE, sucra = "a"), m)
  expect_error(posterior_ranks(smk_fit_RE, sucra = 1), m)
  expect_error(posterior_ranks(smk_fit_RE, sucra = list()), m)
  expect_error(posterior_ranks(smk_fit_RE, sucra = NA), m)
  expect_error(posterior_ranks(smk_fit_RE, sucra = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(posterior_ranks(smk_fit_RE, newdata = "a"), m)
  expect_error(posterior_ranks(smk_fit_RE, newdata = 1), m)
  expect_error(posterior_ranks(smk_fit_RE, newdata = list()), m)
  expect_error(posterior_ranks(smk_fit_RE, newdata = NA), m)
})

test_that("lower_better argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_ranks(smk_fit_RE, lower_better = "a"), m)
  expect_error(posterior_ranks(smk_fit_RE, lower_better = 1), m)
  expect_error(posterior_ranks(smk_fit_RE, lower_better = list()), m)
  expect_error(posterior_ranks(smk_fit_RE, lower_better = NA), m)
  expect_error(posterior_ranks(smk_fit_RE, lower_better = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(posterior_rank_probs(smk_fit_RE, newdata = "a"), m)
  expect_error(posterior_rank_probs(smk_fit_RE, newdata = 1), m)
  expect_error(posterior_rank_probs(smk_fit_RE, newdata = list()), m)
  expect_error(posterior_rank_probs(smk_fit_RE, newdata = NA), m)
})

test_that("lower_better argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_rank_probs(smk_fit_RE, lower_better = "a"), m)
  expect_error(posterior_rank_probs(smk_fit_RE, lower_better = 1), m)
  expect_error(posterior_rank_probs(smk_fit_RE, lower_better = list()), m)
  expect_error(posterior_rank_probs(smk_fit_RE, lower_better = NA), m)
  expect_error(posterior_rank_probs(smk_fit_RE, lower_better = NULL), m)
})

test_that("cumulative argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_rank_probs(smk_fit_RE, cumulative = "a"), m)
  expect_error(posterior_rank_probs(smk_fit_RE, cumulative = 1), m)
  expect_error(posterior_rank_probs(smk_fit_RE, cumulative = list()), m)
  expect_error(posterior_rank_probs(smk_fit_RE, cumulative = NA), m)
  expect_error(posterior_rank_probs(smk_fit_RE, cumulative = NULL), m)
})

test_that("sucra argument", {
  m <- "should be TRUE or FALSE"
  expect_error(posterior_rank_probs(smk_fit_RE, sucra = "a"), m)
  expect_error(posterior_rank_probs(smk_fit_RE, sucra = 1), m)
  expect_error(posterior_rank_probs(smk_fit_RE, sucra = list()), m)
  expect_error(posterior_rank_probs(smk_fit_RE, sucra = NA), m)
  expect_error(posterior_rank_probs(smk_fit_RE, sucra = NULL), m)
})

skip_on_cran()  # Reduce CRAN check time

# Only small number of samples to test
smk_fit_RE <- suppressWarnings(nma(smk_net,
                                   trt_effects = "random",
                                   prior_intercept = normal(scale = 100),
                                   prior_trt = normal(scale = 100),
                                   prior_het = normal(scale = 5),
                                   iter = 10))

test_that(".trt column is correct", {
  rk <- tibble::as_tibble(posterior_ranks(smk_fit_RE))
  expect_identical(paste0("rank[", rk$.trt, "]"),
                   rk$parameter)

  rp <- tibble::as_tibble(posterior_rank_probs(smk_fit_RE))
  expect_identical(paste0("d[", rp$.trt, "]"),
                   rp$parameter)
})


pso_net <- set_ipd(plaque_psoriasis_ipd[complete.cases(plaque_psoriasis_ipd), ],
                   studyc, trtc,
                   r = pasi75)

# Only small number of samples to test
pso_fit <- suppressWarnings(nma(pso_net,
                                trt_effects = "fixed",
                                regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
                                prior_intercept = normal(scale = 10),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                init_r = 0.1,
                                iter = 10))

test_that(".study, .trt columns are correct", {
  rk <- tibble::as_tibble(posterior_ranks(pso_fit))
  expect_identical(paste0("rank[", rk$.study, ": ", rk$.trt, "]"),
                   rk$parameter)

  rp <- tibble::as_tibble(posterior_rank_probs(pso_fit))
  expect_identical(paste0("d[", rp$.study, ": ", rp$.trt, "]"),
                   rp$parameter)
})
