library(multinma)
library(dplyr)

test_that("set_* produces empty nma_data objects", {
  empty_nma_data <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = NULL,
         studies = NULL), class = "nma_data")

  expect_equal(set_ipd(smoking[0, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_arm(smoking[0, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_contrast(smoking[0, ], "studyn", "trtc"), empty_nma_data)
})

test_that("set_* error if data does not inherit data.frame", {
  vec <- 1:5
  msg <- "Argument `data` should be a data frame"

  expect_error(set_ipd(vec), msg)
  expect_error(set_agd_arm(vec), msg)
  expect_error(set_agd_contrast(vec), msg)
})

test_that("set_* error if study not given", {
  expect_error(set_ipd(smoking), "Specify `study`")
  expect_error(set_agd_arm(smoking), "Specify `study`")
  expect_error(set_agd_contrast(smoking), "Specify `study`")
})

test_that("set_* error if trt not given", {
  expect_error(set_ipd(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_arm(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_contrast(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_contrast(smoking, "studyn", "trtc"), "Specify `trt_b`")
})

# Dummy data
agd_arm <- tibble(
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
 disc_neg = -disc,
 bin = rbinom(5, 1, 0.5)
 #Surv =
)

test_that("set_agd_arm - continuous outcome checks work", {
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont), "Specify standard error")
  expect_warning(set_agd_arm(agd_arm, "studyn", "trtc", se = cont_pos), "Ignoring standard error")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_neg), "must be positive")
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_pos)$agd_arm[, c(".y", ".se")],
    transmute(agd_arm, .y = cont, .se = cont_pos))
})

test_that("set_agd_arm - count outcome checks work", {
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
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = cont_pos)$agd_arm[, c(".r", ".n", ".E")],
    transmute(agd_arm, .r = disc, .n = disc_p1, .E = cont_pos))
})

test_that("set_ipd - continuous outcome checks work", {
  expect_error(set_ipd(agd_arm, "studyn", "trtc", y = trtc), "must be numeric")
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", y = cont)$ipd[, ".y"],
    transmute(agd_arm, .y = cont))
})

test_that("set_ipd - binary outcome checks work", {
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = trtc), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = cont), "must equal 0 or 1")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = disc_neg), "must equal 0 or 1")
  expect_warning(set_ipd(agd_arm, "studyn", "trtc", E = cont_pos), "Ignoring `E`")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cont_neg), "must be positive")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = trtc), "must be numeric")
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cont_pos)$ipd[, c(".r", ".E")],
    transmute(agd_arm, .r = bin, .E = cont_pos))
})

# Dummy contrast data
agd_contrast <- agd_arm %>%
  group_by(studyn) %>%
  mutate(trt_b = first(trtf)) %>%
  filter(trt_b != trtf) %>%
  ungroup()

test_that("set_agd_contrast - continuous outcome checks work", {
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc"), "Specify `trt_b`")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont), "Specify standard error")
  expect_warning(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", se = cont_pos), "Ignoring standard error")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = cont_neg), "must be positive")
  expect_equivalent(
    set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = cont_pos)$agd_contrast[, c(".y", ".se")],
    transmute(agd_contrast, .y = cont, .se = cont_pos))
})
