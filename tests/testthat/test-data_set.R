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

test_that("set_* error if study not given or missing values", {
  expect_error(set_ipd(smoking), "Specify `study`")
  expect_error(set_agd_arm(smoking), "Specify `study`")
  expect_error(set_agd_contrast(smoking), "Specify `study`")

  smk_miss <- smoking
  smk_miss[1, "studyn"] <- NA
  expect_error(set_ipd(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_arm(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_contrast(smk_miss, "studyn", "trtc", "trtn"), "cannot contain missing values")
})

test_that("set_* error if trt not given or missing values", {
  expect_error(set_ipd(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_arm(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_contrast(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_contrast(smoking, "studyn", "trtc"), "Specify `trt_b`")

  smk_miss <- smoking
  smk_miss[1, "trtc"] <- NA
  expect_error(set_ipd(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_arm(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_contrast(smk_miss, "studyn", "trtn", "trtc"), "cannot contain missing values")
})

# Dummy data
agd_arm <- tibble(
  studyn = c(1, 1, 2, 2, 2),
  studyc = letters[studyn],
  studyf = factor(studyc),
  trtn = c(1, 2, 1, 2, 3),
  trtc = LETTERS[trtn],
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
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", se = cont_pos), "Specify continuous outcome `y`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_neg), "must be positive")
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_pos)$agd_arm[, c(".y", ".se")],
    transmute(agd_arm, .y = cont, .se = cont_pos))
})

test_that("set_agd_arm - count outcome checks work", {
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc), "Specify denominator")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", n = disc_p1), "Specify outcome count `r`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = trtc, n = disc_p1), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = cont, n = disc), "must be integer")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = cont), "must be integer")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc_neg, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc_p1, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_neg), "greater than zero")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", E = cont_pos), "Specify outcome count `r`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = cont_neg), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = trtc), "must be numeric")
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1)$agd_arm[, c(".r", ".n")],
    transmute(agd_arm, .r = disc, .n = disc_p1))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, E = cont_pos)$agd_arm[, c(".r", ".E")],
    transmute(agd_arm, .r = disc, .E = cont_pos))
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
  expect_error(set_ipd(agd_arm, "studyn", "trtc", E = cont_pos), "Specify count `r`")
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
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", se = cont_pos), "Specify continuous outcome `y`")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = cont_neg), "must be positive")
  expect_equivalent(
    set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b", y = cont, se = cont_pos)$agd_contrast[, c(".y", ".se")],
    transmute(agd_contrast, .y = cont, .se = cont_pos))
})

test_that("set_* - take one and only one outcome", {
  m <- "specify one and only one outcome"
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, y = cont), m)
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, y = cont, se = cont_pos), m)
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", "trt_b"), m)
})

test_that("set_* `.trt` column is correct", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, trt_b, y = cont, se = cont_pos)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, trt_b, y = cont, se = cont_pos)$agd_contrast$.trt_b,
               agd_contrast$trt_b)
})

test_that("set_* `.study` column is correct", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, trt_b, y = cont, se = cont_pos)$agd_contrast$.study,
               agd_contrast$studyf)
})

test_that("set_* return `treatments` factor", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$treatments,
               factor(LETTERS[1:3]))
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$treatments,
               factor(LETTERS[1:3]))
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, trt_b, y = cont, se = cont_pos)$treatments,
               factor(LETTERS[1:3]))
})

test_that("set_* return `studies` factor", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$studies,
               factor(letters[1:2]))
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$studies,
               factor(letters[1:2]))
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, trt_b, y = cont, se = cont_pos)$studies,
               factor(letters[1:2]))
})

make_na <- function(x, n) {
  x[sample.int(length(x), n)] <- NA
  return(x)
}

test_that("set_* error if outcomes contain missing values", {
  agd_arm_miss <- agd_arm %>% mutate_at(vars(cont:bin), ~make_na(., 1))
  agd_contrast_miss <- agd_contrast %>% mutate_at(vars(cont:bin), ~make_na(., 1))

  m <- "contains missing values"

  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, y = cont, se = cont_pos), m)
  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, r = disc, n = disc_p1), m)
  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, r = disc, E = cont_pos), m)

  expect_error(set_ipd(agd_arm_miss, studyn, trtn, y = cont), m)
  expect_error(set_ipd(agd_arm_miss, studyn, trtn, r = disc), m)
  expect_error(set_ipd(agd_arm_miss, studyn, trtn, r = disc, E = cont_pos), m)

  expect_error(set_agd_contrast(agd_contrast_miss, studyn, trtn, trt_b, y = cont, se = cont_pos), m)
})
