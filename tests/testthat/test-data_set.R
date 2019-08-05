library(multinma)

test_that("set_* produces empty nma_data objects", {
  empty_nma_data <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = NULL,
         studies = NULL), class = "nma_data")

  expect_equal(set_ipd(smoking[NA, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_arm(smoking[NA, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_contrast(smoking[NA, ], "studyn", "trtc"), empty_nma_data)
})

test_that("set_* error if data does not inherit data.frame", {
  vec <- 1:5
  msg <- "Argument `data` should be a data frame"

  expect_error(set_ipd(vec), msg)
  expect_error(set_agd_arm(vec), msg)
  expect_error(set_agd_contrast(vec), msg)
})

# Dummy data
agd_arm <- tibble::tibble(
 studyn = c(1, 1, 2, 2, 2),
 studyc = c("a", "a", "b", "b", "b"),
 studyf = factor(studyc),
 trtn = c(1, 2, 1, 2, 3),
 trtc = c("A", "B", "A", "B", "C"),
 trtf = factor(trtc),
 cont = rnorm(5),
 cont_pos = abs(cont),
 cont_neg = -cont_pos,
 disc = rbinom(5, 20, 0.5) + 1,
 disc_p1 = disc + 1,
 disc_m1 = disc - 1,
 disc_neg = -disc
 #Surv =
)

test_that("continuous outcome checks work", {
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont), "Specify standard error")
  expect_warning(set_agd_arm(agd_arm, "studyn", "trtc", se = cont_pos), "Ignoring standard error")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_neg), "must be positive")
})

test_that("count outcome checks work", {
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc), "Specify denominator")
  expect_warning(set_agd_arm(agd_arm, "studyn", "trtc", n = disc_p1), "Ignoring `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = trtc, n = disc_p1), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = cont, n = disc), "must be integer")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = cont), "must be integer")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc_neg, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc_p1, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_neg), "greater than zero")
  expect_warning(set_agd_arm(agd_arm, "studyn", "trtc", E = cont_pos), "Ignoring `E`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = cont_neg), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = trtc), "must be numeric")
})
