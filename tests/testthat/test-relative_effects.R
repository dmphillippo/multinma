
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

test_that("all_contrasts argument", {
  m <- "should be TRUE or FALSE"
  expect_error(relative_effects(smk_fit_RE, all_contrasts = "a"), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = 1), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = list()), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = NA), m)
  expect_error(relative_effects(smk_fit_RE, all_contrasts = NULL), m)
})

test_that("trt_ref argument", {
  m <- "does not match a treatment in the network"
  expect_error(relative_effects(smk_fit_RE, trt_ref = "a"), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = 1), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = list("a")), m)
  expect_error(relative_effects(smk_fit_RE, trt_ref = NA), m)

  # expect_warning(relative_effects(smk_fit_RE, all_contrasts = TRUE, trt_ref = "Self-help"), "Ignoring `trt_ref`")
})

test_that("probs argument", {
  m <- "numeric vector of probabilities"
  expect_error(relative_effects(smk_fit_RE, probs = "a"), m)
  expect_error(relative_effects(smk_fit_RE, probs = -1), m)
  expect_error(relative_effects(smk_fit_RE, probs = 1.5), m)
  expect_error(relative_effects(smk_fit_RE, probs = Inf), m)
  expect_error(relative_effects(smk_fit_RE, probs = list()), m)
  expect_error(relative_effects(smk_fit_RE, probs = NA), m)
  expect_error(relative_effects(smk_fit_RE, probs = NULL), m)
})

test_that("summary argument", {
  m <- "should be TRUE or FALSE"
  expect_error(relative_effects(smk_fit_RE, summary = "a"), m)
  expect_error(relative_effects(smk_fit_RE, summary = 1), m)
  expect_error(relative_effects(smk_fit_RE, summary = list()), m)
  expect_error(relative_effects(smk_fit_RE, summary = NA), m)
  expect_error(relative_effects(smk_fit_RE, summary = NULL), m)
})

test_that("newdata argument", {
  m <- "not a data frame"
  expect_error(relative_effects(smk_fit_RE, newdata = "a"), m)
  expect_error(relative_effects(smk_fit_RE, newdata = 1), m)
  expect_error(relative_effects(smk_fit_RE, newdata = list()), m)
  expect_error(relative_effects(smk_fit_RE, newdata = NA), m)
})


skip_on_cran()  # Reduce CRAN check time

# Only small number of samples to test
smk_fit_RE <- suppressWarnings(nma(smk_net,
                                   trt_effects = "random",
                                   prior_intercept = normal(scale = 100),
                                   prior_trt = normal(scale = 100),
                                   prior_het = normal(scale = 5),
                                   iter = 10))

test_that(".trta, .trtb columns are correct", {
  re1 <- tibble::as_tibble(relative_effects(smk_fit_RE))
  expect_identical(paste0("d[", re1$.trtb, "]"),
                   re1$parameter)
  expect_identical(as.character(re1$.trta), rep("No intervention", times = nrow(re1)))

  re2 <- tibble::as_tibble(relative_effects(smk_fit_RE, all_contrasts = TRUE))
  expect_identical(paste0("d[", re2$.trtb, " vs. ", re2$.trta, "]"),
                   re2$parameter)

  re3 <- tibble::as_tibble(relative_effects(smk_fit_RE, trt_ref = "Self-help"))
  expect_identical(paste0("d[", re3$.trtb, "]"),
                   re3$parameter)
  expect_identical(as.character(re3$.trta), rep("Self-help", times = nrow(re3)))
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

pso_new <- data.frame(durnpso = 10, prevsys = TRUE, bsa = 20, weight = 75, psa = FALSE, study = c("One", "Two"))

test_that("newdata validation", {
  expect_error(relative_effects(pso_fit, newdata = pso_new[-1], study = study),
               'Regression variable "durnpso" not found in `newdata`')
  expect_error(relative_effects(pso_fit, newdata = pso_new[-(1:2)], study = study),
               'Regression variables "durnpso" and "prevsys" not found in `newdata`')

  make_bad <- function(x, vars = "durnpso") {
    bad <- pso_new
    bad[vars] <- x
    return(bad)
  }

  expect_error(relative_effects(pso_fit, newdata = make_bad(NA), study = study),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(relative_effects(pso_fit, newdata = make_bad(NaN), study = study),
               "missing or infinite values in `newdata`: durnpso")
  expect_error(relative_effects(pso_fit, newdata = make_bad(Inf), study = study),
               "missing or infinite values in `newdata`: durnpso")

  expect_error(relative_effects(pso_fit, newdata = make_bad(NA, c("durnpso", "psa")), study = study),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(relative_effects(pso_fit, newdata = make_bad(NaN, c("durnpso", "psa")), study = study),
               "missing or infinite values in `newdata`: durnpso, psa")
  expect_error(relative_effects(pso_fit, newdata = make_bad(Inf, c("durnpso", "psa")), study = study),
               "missing or infinite values in `newdata`: durnpso, psa")
})

test_that(".study, .trta, .trtb columns are correct", {
  re1 <- tibble::as_tibble(relative_effects(pso_fit))
  expect_identical(paste0("d[", re1$.study, ": ", re1$.trtb, "]"),
                   re1$parameter)
  expect_identical(as.character(re1$.trta), rep(levels(pso_net$treatments)[1], times = nrow(re1)))

  re2 <- tibble::as_tibble(relative_effects(pso_fit, all_contrasts = TRUE))
  expect_identical(paste0("d[", re2$.study, ": ", re2$.trtb, " vs. ", re2$.trta, "]"),
                   re2$parameter)

  re3 <- tibble::as_tibble(relative_effects(pso_fit, trt_ref = "ETN"))
  expect_identical(paste0("d[", re3$.study, ": ", re3$.trtb, "]"),
                   re3$parameter)
  expect_identical(as.character(re3$.trta), rep("ETN", times = nrow(re3)))

  re4 <- tibble::as_tibble(relative_effects(pso_fit, newdata = pso_new, study = study))
  expect_identical(paste0("d[", re4$.study, ": ", re4$.trtb, "]"),
                   re4$parameter)
  expect_identical(as.character(re4$.trta), rep(levels(pso_net$treatments)[1], times = nrow(re4)))

  re5 <- tibble::as_tibble(relative_effects(pso_fit, all_contrasts = TRUE, newdata = pso_new, study = study))
  expect_identical(paste0("d[", re5$.study, ": ", re5$.trtb, " vs. ", re5$.trta, "]"),
                   re5$parameter)

  re6 <- tibble::as_tibble(relative_effects(pso_fit, trt_ref = "ETN", newdata = pso_new, study = study))
  expect_identical(paste0("d[", re6$.study, ": ", re6$.trtb, "]"),
                   re6$parameter)
  expect_identical(as.character(re6$.trta), rep("ETN", times = nrow(re6)))
})
