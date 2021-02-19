
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

test_that("baseline argument", {
  m <- "should be specified using distr"
  expect_error(predict(smk_fit_RE, baseline = "a"), m)
  expect_error(predict(smk_fit_RE, baseline = 1), m)
  expect_error(predict(smk_fit_RE, baseline = list("a")), m)
  expect_error(predict(smk_fit_RE, baseline = NA), m)
})

test_that("trt_ref argument", {
  m <- "does not match a treatment in the network"
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = "a"), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = 1), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = list("a")), m)
  expect_error(predict(smk_fit_RE, baseline = distr(qnorm, mean = -1, sd = 0.01), trt_ref = NA), m)
})

test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(predict(smk_fit_RE, probs = "a"), m)
  expect_error(predict(smk_fit_RE, probs = -1), m)
  expect_error(predict(smk_fit_RE, probs = 1.5), m)
  expect_error(predict(smk_fit_RE, probs = Inf), m)
  expect_error(predict(smk_fit_RE, probs = list()), m)
  expect_error(predict(smk_fit_RE, probs = NA), m)
  expect_error(predict(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(predict(smk_fit_RE, summary = "a"), m)
  expect_error(predict(smk_fit_RE, summary = 1), m)
  expect_error(predict(smk_fit_RE, summary = list()), m)
  expect_error(predict(smk_fit_RE, summary = NA), m)
  expect_error(predict(smk_fit_RE, summary = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(predict(smk_fit_RE, newdata = "a"), m)
  expect_error(predict(smk_fit_RE, newdata = 1), m)
  expect_error(predict(smk_fit_RE, newdata = list()), m)
  expect_error(predict(smk_fit_RE, newdata = NA), m)
})

test_that("type argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, type = "a"), m)
  expect_error(predict(smk_fit_RE, type = "lin"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, type = 1), m2)
  expect_error(predict(smk_fit_RE, type = list("a")), m2)
  expect_error(predict(smk_fit_RE, type = NA), m2)
})

test_that("level argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, level = "a"), m)
  expect_error(predict(smk_fit_RE, level = "agg"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, level = 1), m2)
  expect_error(predict(smk_fit_RE, level = list("a")), m2)
  expect_error(predict(smk_fit_RE, level = NA), m2)

  expect_error(predict(smk_fit_RE, level = "individual"),
               "Cannot produce individual predictions without a regression model.")
})

test_that("baseline_type argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, baseline_type = "a"), m)
  expect_error(predict(smk_fit_RE, baseline_type = "lin"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, baseline_type = 1), m2)
  expect_error(predict(smk_fit_RE, baseline_type = list("a")), m2)
  expect_error(predict(smk_fit_RE, baseline_type = NA), m2)
})

test_that("baseline_level argument", {
  m <- "must be one of"
  expect_error(predict(smk_fit_RE, baseline_level = "a"), m)
  expect_error(predict(smk_fit_RE, baseline_level = "agg"), m)

  m2 <- "must be a character vector"
  expect_error(predict(smk_fit_RE, baseline_level = 1), m2)
  expect_error(predict(smk_fit_RE, baseline_level = list("a")), m2)
  expect_error(predict(smk_fit_RE, baseline_level = NA), m2)
})


skip_on_cran()  # Reduce CRAN check time

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

pso_new <- data.frame(durnpso = 10, prevsys = TRUE, bsa = 20, weight = 75, psa = FALSE, study = c("One", "Two"))

test_that("baseline and newdata for regression models", {
  m <- "Specify both `newdata` and `baseline`, or neither"
  expect_error(predict(pso_fit, newdata = pso_new), m)
  expect_error(predict(pso_fit, baseline = distr(qnorm, 1, 1)), m)
})

test_that("baseline argument", {
  m <- "should be a single distr\\(\\) specification, a list of distr\\(\\) specifications, or NULL"
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = "a"), m)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = 1), m)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list("a")), m)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = NA), m)

  m2 <- "or a list of length 2"
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(distr(qnorm, 1, 1),
                                                                                  distr(qnorm, 2, 1),
                                                                                  distr(qnorm, 3, 1))), m2)
  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(One = distr(qnorm, 1, 1),
                                                                                  Two = distr(qnorm, 2, 1),
                                                                                  Three = distr(qnorm, 3, 1))), m2)

  expect_error(predict(pso_fit, study = study, newdata = pso_new, baseline = list(One = distr(qnorm, 1, 1),
                                                                                  Three = distr(qnorm, 3, 1))),
               "must match all study names")
})

test_that("newdata validation", {
  expect_error(predict(pso_fit, newdata = pso_new[-1], study = study, baseline = distr(qnorm, 1, 1)),
               'Regression variable "durnpso" not found in `newdata`')
  expect_error(predict(pso_fit, newdata = pso_new[-(1:2)], study = study, baseline = distr(qnorm, 1, 1)),
               'Regression variables "durnpso" and "prevsys" not found in `newdata`')

  make_bad <- function(x, vars = "durnpso") {
    bad <- pso_new
    bad[vars] <- x
    return(bad)
  }

  expect_error(predict(pso_fit, newdata = make_bad(NA), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(predict(pso_fit, newdata = make_bad(NaN), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(predict(pso_fit, newdata = make_bad(Inf), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso")

  expect_error(predict(pso_fit, newdata = make_bad(NA, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(predict(pso_fit, newdata = make_bad(NaN, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(predict(pso_fit, newdata = make_bad(Inf, c("durnpso", "psa")), study = study, baseline = distr(qnorm, 1, 1)),
               "missing or infinite values in `newdata`: durnpso, psa")
})
