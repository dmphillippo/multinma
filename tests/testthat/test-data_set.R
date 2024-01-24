library(multinma)
library(dplyr)

test_that("set_* produces empty nma_data objects", {
  empty_nma_data <- structure(
    list(agd_arm = NULL,
         agd_contrast = NULL,
         ipd = NULL,
         treatments = NULL,
         classes = NULL,
         studies = NULL), class = "nma_data")

  expect_equal(set_ipd(smoking[0, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_arm(smoking[0, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_contrast(smoking[0, ], "studyn", "trtc"), empty_nma_data)
  expect_equal(set_agd_surv(smoking[0, ], "studyn", "trtc"), empty_nma_data)
})

test_that("set_* error if data does not inherit data.frame", {
  vec <- 1:5
  msg <- "Argument `data` should be a data frame"

  expect_error(set_ipd(vec), msg)
  expect_error(set_agd_arm(vec), msg)
  expect_error(set_agd_contrast(vec), msg)
  expect_error(set_agd_surv(vec), msg)
})

test_that("set_* error if study not given, missing values, or not regular 1D column", {
  expect_error(set_ipd(smoking), "Specify `study`")
  expect_error(set_agd_arm(smoking), "Specify `study`")
  expect_error(set_agd_contrast(smoking), "Specify `study`")
  expect_error(set_agd_surv(smoking), "Specify `study`")

  smk_miss <- smoking
  smk_miss[1, "studyn"] <- NA
  expect_error(set_ipd(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_arm(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_contrast(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_surv(smk_miss, "studyn", "trtc"), "cannot contain missing values")

  expect_error(set_ipd(smoking, cbind(studyn, studyn)), "must be a regular column")
  expect_error(set_ipd(smoking, list(studyn)), "must be a regular column")
  expect_error(set_agd_arm(smoking, cbind(studyn, studyn)), "must be a regular column")
  expect_error(set_agd_arm(smoking, list(studyn)), "must be a regular column")
  expect_error(set_agd_contrast(smoking, cbind(studyn, studyn)), "must be a regular column")
  expect_error(set_agd_contrast(smoking, list(studyn)), "must be a regular column")
  expect_error(set_agd_surv(smoking, cbind(studyn, studyn)), "must be a regular column")
  expect_error(set_agd_surv(smoking, list(studyn)), "must be a regular column")
})

test_that("set_* error if single-arm studies included", {
  m <- "Single-arm studies are not supported"

  s <- tibble(study = c("a", "b", "b", "c"), trt = c("A", "A", "A", "B"), r = 1, n = 2, time = 1, status = 1)
  expect_error(set_ipd(s, study, trt, r = r), paste0(m, '.+studies "a", "b" and "c"'))
  expect_error(set_agd_arm(s, study, trt, r = r, n = n), paste0(m, '.+studies "a" and "c"'))
  expect_error(set_agd_contrast(s, study, trt, y = r, se = n), paste0(m, '.+studies "a" and "c"'))

  # Allowed with message for survival outcomes
  expect_message(set_ipd(s, study, trt, Surv = Surv(time, status)), 'Single-arm studies present in the network: "a", "b" and "c"')
  expect_message(set_agd_surv(s, study, trt, Surv = Surv(time, status)), 'Single-arm studies present in the network: "a", "b" and "c"')
})

test_that("set_* error if trt not given, missing values, or not regular 1D column", {
  expect_error(set_ipd(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_arm(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_contrast(smoking, "studyn"), "Specify `trt`")
  expect_error(set_agd_surv(smoking, "studyn"), "Specify `trt`")

  smk_miss <- smoking
  smk_miss[1, "trtc"] <- NA
  expect_error(set_ipd(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_arm(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_contrast(smk_miss, "studyn", "trtc"), "cannot contain missing values")
  expect_error(set_agd_surv(smk_miss, "studyn", "trtc"), "cannot contain missing values")

  expect_error(set_ipd(smoking, studyn, cbind(trtc, trtc)), "must be a regular column")
  expect_error(set_ipd(smoking, studyn, list(trtc)), "must be a regular column")
  expect_error(set_agd_arm(smoking, studyn, cbind(trtc, trtc)), "must be a regular column")
  expect_error(set_agd_arm(smoking, studyn, list(trtc)), "must be a regular column")
  expect_error(set_agd_contrast(smoking, studyn, cbind(trtc, trtc)), "must be a regular column")
  expect_error(set_agd_contrast(smoking, studyn, list(trtc)), "must be a regular column")
  expect_error(set_agd_surv(smoking, studyn, cbind(trtc, trtc)), "must be a regular column")
  expect_error(set_agd_surv(smoking, studyn, list(trtc)), "must be a regular column")
})

# Dummy data
agd_arm <- tibble(
  studyn = c(1, 1, 2, 2, 2),
  studyc = letters[studyn],
  studyf = factor(studyc),
  studyf2 = forcats::fct_rev(studyf),
  trtn = c(1, 2, 1, 2, 3),
  trtc = LETTERS[trtn],
  trtf = factor(trtc),
  trtf2 = forcats::fct_rev(trtf),
  tclassn = c(1, 2, 1, 2, 2),
  tclassc = letters[tclassn],
  tclassf = factor(tclassc),
  tclassf2 = forcats::fct_rev(tclassf),
  cont = rnorm(5),
  cont_pos = abs(cont),
  cont_neg = -cont_pos,
  cont_inf = c(cont_pos[1:4], Inf),
  cont_nan = c(cont_pos[1:4], NaN),
  disc = rbinom(5, 20, 0.5) + 1,
  disc_p1 = disc + 1,
  disc_m1 = disc - 1,
  disc_neg = -disc,
  disc_inf = c(disc[1:4], Inf),
  disc_nan = c(disc[1:4], NaN),
  disc_na = c(disc[1:4], NA),
  bin = sample(c(0,0,0,0,1,1,1,1), 5)
  #Surv =
)

test_that("set_agd_arm - continuous outcome checks work", {
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont), "Specify standard error")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", se = cont_pos), "Specify continuous outcome `y`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = trtc, se = cont_pos), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = as.character(trtn), se = cont_pos), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = as.character(trtn)), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_neg), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = -cont_pos), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_inf), "cannot be infinite")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", y = cont, se = cont_nan), "cannot be NaN")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cbind(cont, cont), se = cont_pos),
               "must be a regular column")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cbind(cont_pos, cont_pos)),
               "must be a regular column")
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
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = -disc, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc_p1, n = disc), "must be between 0 and `n`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_neg), "greater than zero")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = -disc), "greater than zero")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", E = cont_pos), "Specify outcome count `r`")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = cont_neg), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = -cont_pos), "must be positive")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = trtc), "must be numeric")
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, E = as.character(trtn)), "must be numeric")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = cbind(disc, disc), n = disc_p1), "must be a regular column")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = disc, n = cbind(disc_p1, disc_p1)), "must be a regular column")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = disc, E = cbind(cont_pos, cont_pos)), "must be a regular column")
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1)$agd_arm[, c(".r", ".n")],
    transmute(agd_arm, .r = disc, .n = disc_p1))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = "disc", n = "disc_p1")$agd_arm[, c(".r", ".n")],
    transmute(agd_arm, .r = disc, .n = disc_p1))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc + 1)$agd_arm[, c(".r", ".n")],
    transmute(agd_arm, .r = disc, .n = disc_p1))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = floor(disc/2), n = disc)$agd_arm[, c(".r", ".n")],
    transmute(agd_arm, .r = floor(disc/2), .n = disc))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = disc, E = cont_pos)$agd_arm[, c(".r", ".E")],
    transmute(agd_arm, .r = disc, .E = cont_pos))
  expect_equivalent(
    set_agd_arm(agd_arm, "studyn", "trtc", r = "disc", E = "cont_pos")$agd_arm[, c(".r", ".E")],
    transmute(agd_arm, .r = disc, .E = cont_pos))
})


multi_inclusive <- tribble(~r_a, ~r_b, ~r_c,
                           1, 1, 1,
                           5, 4, 1,
                           5, 2, 2,
                           10, 5, 0,
                           6, 0, 0)
agd_arm_multi_i <- bind_cols(agd_arm, multi_inclusive)

multi_exclusive <- tribble(~r_a, ~r_b, ~r_c,
                           0, 0, 1,
                           1, 3, 1,
                           3, 0, 2,
                           5, 5, 0,
                           6, 0, 0)
agd_arm_multi_e <- bind_cols(agd_arm, multi_exclusive)

test_that("set_agd_arm - multinomial outcome checks work", {
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc)),
               "At least 2 outcomes", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, disc)),
               "Duplicate outcome category labels", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, 1:2)),
               "must be the same length", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, disc_inf)),
               "cannot be Inf", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, disc_nan)),
               "cannot be NaN", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, trtc)),
               "must be numeric", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, cont_pos)),
               "must be integer", class = "error")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, r = multi(disc, disc_neg)),
               "must be non-negative", class = "error")

  expect_error(set_agd_arm(agd_arm_multi_i, studyn, trtc, r = multi(r_c, r_b, r_a, inclusive = TRUE)),
               "must be decreasing or constant", class = "error")
  expect_equivalent(unclass(set_agd_arm(agd_arm_multi_i, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = TRUE))$agd_arm$.r),
                    as.matrix(multi_exclusive))

  expect_equivalent(unclass(set_agd_arm(agd_arm_multi_e, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = FALSE))$agd_arm$.r),
                    as.matrix(multi_exclusive))
})

test_that("set_agd_arm - sample size checks work", {
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = trtc),
               "must be numeric")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = as.character(trtn)),
               "must be numeric")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = cont),
               "must be integer")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = disc_neg),
               "must be greater than zero")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = -disc),
               "must be greater than zero")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = disc_inf),
               "cannot be infinite")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = disc_nan),
               "cannot be NaN")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = cbind(disc, disc)),
               "must be a regular column")
  expect_error(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos, sample_size = list(disc)),
               "must be a regular column")
  expect_message(set_agd_arm(agd_arm, studyn, trtc, y = cont, se = cont_pos), "`sample_size` not provided")

  expect_equal(set_agd_arm(agd_arm, studyn, trtc, r = disc, n = disc_p1)$agd_arm$.sample_size, agd_arm$disc_p1)
  expect_equal(set_agd_arm(agd_arm, studyn, trtc, r = floor(disc/2), n = disc, sample_size = disc + 1)$agd_arm$.sample_size, agd_arm$disc_p1)
  expect_equal(set_agd_arm(agd_arm_multi_i, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = TRUE))$agd_arm$.sample_size, multi_inclusive$r_a)
  expect_equal(set_agd_arm(agd_arm_multi_e, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = FALSE))$agd_arm$.sample_size, multi_inclusive$r_a)
})

test_that("set_ipd - continuous outcome checks work", {
  expect_error(set_ipd(agd_arm, "studyn", "trtc", y = trtc), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", y = as.character(trtn)), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", y = cbind(cont, cont)), "must be a regular column")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", y = list(cont)), "must be a regular column")
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", y = cont)$ipd[, ".y"],
    transmute(agd_arm, .y = cont))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", y = "cont")$ipd[, ".y"],
    transmute(agd_arm, .y = cont))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", y = cont/2)$ipd[, ".y"],
    transmute(agd_arm, .y = cont/2))
})

test_that("set_ipd - binary outcome checks work", {
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = trtc), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = as.character(trtn)), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = cont), "must equal 0 or 1")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = disc_neg), "must equal 0 or 1")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin + 1), "must equal 0 or 1")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", E = cont_pos), "Specify count `r`")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cont_neg), "must be positive")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = -cont_pos), "must be positive")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = trtc), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = as.character(trtn)), "must be numeric")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = cbind(bin, bin)), "must be a regular column")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = cbind(bin, bin), E = cont_pos), "must be a regular column")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cbind(cont_pos, cont_pos)), "must be a regular column")
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cont_pos)$ipd[, c(".r", ".E")],
    transmute(agd_arm, .r = bin, .E = cont_pos))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", r = "bin", E = "cont_pos")$ipd[, c(".r", ".E")],
    transmute(agd_arm, .r = bin, .E = cont_pos))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", r = bin, E = cont_pos/2)$ipd[, c(".r", ".E")],
    transmute(agd_arm, .r = bin, .E = cont_pos/2))
})

test_that("set_ipd - multinomial outcome checks work", {
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(bin)),
               "At least 2 outcomes", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(bin, bin)),
               "Duplicate outcome category labels", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(bin, 1:2)),
               "must be the same length", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(disc, disc_inf)),
               "cannot be Inf", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(disc, disc_nan)),
               "cannot be NaN", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(disc, trtc)),
               "must be numeric", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(bin, bin + 1)),
               "must equal 0 or 1", class = "error")
  expect_error(set_ipd(agd_arm, studyn, trtc, r = multi(bin, -bin)),
               "must be non-negative", class = "error")

  i_multi_inclusive <- tribble(~r_a, ~r_b, ~r_c,
                               1, 1, 1,
                               1, 1, 0,
                               1, 0, 0,
                               1, 1, NA,
                               1, NA, 0)
  ipd_multi_i <- bind_cols(agd_arm, i_multi_inclusive)

  i_multi_exclusive <- tribble(~r_a, ~r_b, ~r_c,
                               0, 0, 1,
                               0, 1, 0,
                               1, 0, 0,
                               0, 1, NA,
                               1, NA, 0)
  ipd_multi_e <- bind_cols(agd_arm, i_multi_exclusive)

  expect_error(set_ipd(ipd_multi_i[-4,], study = I(1), trtc, r = multi(r_c, r_b, r_a, inclusive = TRUE)),
               "must be decreasing or constant", class = "error")
  expect_error(set_ipd(ipd_multi_i[-c(4,5),], study = I(1), trtc, r = multi(r_b, r_c, inclusive = TRUE)),
               "Individual without outcomes in any category, row 3", class = "error")
  expect_equivalent(unclass(set_ipd(ipd_multi_i, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = TRUE))$ipd$.r),
                    as.matrix(i_multi_exclusive))

  expect_error(set_ipd(ipd_multi_e[-c(4,5),], study = I(1), trtc, r = multi(r_b, r_c, inclusive = FALSE)),
               "Individual without outcomes in any category, row 3", class = "error")
  expect_error(set_ipd(ipd_multi_e, studyn, trtc, r = multi(c = 1, r_a, r_b, inclusive = FALSE)),
               "Individuals with outcomes in more than one category, rows 2, 3, 4 and 5", class = "error")
  expect_equivalent(unclass(set_ipd(ipd_multi_e, studyn, trtc, r = multi(r_a, r_b, r_c, inclusive = FALSE))$ipd$.r),
                    as.matrix(i_multi_exclusive))
})

test_that("set_ipd - survival outcome checks work", {
  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = trtc), "must be a `Surv` object")
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, bin))$ipd[[".Surv"]],
    with(agd_arm, Surv(cont_pos, bin)))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos + 1, bin))$ipd[[".Surv"]],
    with(agd_arm, Surv(cont_pos, cont_pos + 1, bin)))
  expect_equivalent(
    set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos + 1, c(0, 1, 3, 3, 3), type = "interval"))$ipd[[".Surv"]],
    with(agd_arm, Surv(cont_pos, cont_pos + 1, c(0, 1, 3, 3, 3), type = "interval")))

  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, bin, type = "mstate")), 'type "mright" is not supported')

  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, bin)), "must have strictly positive outcome times")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_inf, bin)), "infinite times")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_nan, bin)), "missing times")

  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, cont_pos, bin)), "must have non-negative start times")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, cont_pos, c(3, 3, 3, 3, 3), type = "interval")), "must have non-negative start times")

  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, rep_len(c(0, 1, NA), nrow(agd_arm)))), "missing event status values")
  expect_error(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv("a", bin)), "not numeric")
  expect_error(suppressWarnings(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos/2, bin))), "missing times")
  expect_error(suppressWarnings(set_ipd(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos*2, bin*5, type = "interval"))), "missing event status values")
})

test_that("set_agd_surv - survival outcome checks work", {
  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = trtc), "must be a `Surv` object")
  expect_equivalent(
    tidyr::unnest(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, bin))$agd_arm, ".Surv")$.Surv,
    with(agd_arm, Surv(cont_pos, bin)))

  expect_equivalent(
    tidyr::unnest(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos + 1, bin))$agd_arm, ".Surv")$.Surv,
    with(agd_arm, Surv(cont_pos, cont_pos + 1, bin)))
  expect_equivalent(
    tidyr::unnest(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos + 1, c(0, 1, 3, 3, 3), type = "interval"))$agd_arm, ".Surv")$.Surv,
    with(agd_arm, Surv(cont_pos, cont_pos + 1, c(0, 1, 3, 3, 3), type = "interval")))

  expect_equivalent(
    tidyr::unnest(set_agd_surv(ndmm_agd, studyf, trtf, Surv = Surv(eventtime, status))$agd_arm, ".Surv")$.Surv,
    with(ndmm_agd, Surv(eventtime, status)))
  expect_equivalent(
    tidyr::unnest(set_agd_surv(ndmm_agd, studyf, trtf,
                               Surv = Surv(eventtime, status),
                               covariates = ndmm_agd_covs)$agd_arm, ".Surv")$.Surv,
    with(ndmm_agd, Surv(eventtime, status)))

  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, bin, type = "mstate")), 'type "mright" is not supported')

  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, bin)), "must have strictly positive outcome times")
  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_inf, bin)), "infinite times")
  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_nan, bin)), "missing times")

  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, cont_pos, bin)), "must have non-negative start times")
  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_neg, cont_pos, c(3, 3, 3, 3, 3), type = "interval")), "must have non-negative start times")

  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, rep_len(c(0, 1, NA), nrow(agd_arm)))), "missing event status values")
  expect_error(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv("a", bin)), "not numeric")
  expect_error(suppressWarnings(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos/2, bin))), "missing times")
  expect_error(suppressWarnings(set_agd_surv(agd_arm, "studyn", "trtc", Surv = Surv(cont_pos, cont_pos*2, bin*5, type = "interval"))), "missing event status values")
})

test_that("set_agd_surv - covariate checks work", {
  expect_error(set_agd_surv(ndmm_agd, studyf, trtf, Surv = Surv(eventtime, status), covariates = list()), "should be a data frame")

  expect_equivalent(
    dplyr::select(set_agd_surv(ndmm_agd, studyf, trtf, Surv = Surv(eventtime, status),
                                  covariates = ndmm_agd_covs)$agd_arm,
                     study, trt, sample_size = .sample_size,
                     age_min:male) %>%
      dplyr::arrange(study, trt),
    dplyr::select(ndmm_agd_covs, -studyf, -trtf))

  expect_error(set_agd_surv(ndmm_agd, studyf, trtf, Surv = Surv(eventtime, status),
                            covariates = ndmm_agd_covs[-1,]),
               "Not all study arms in `data` have matching rows in `covariates`")
})

# Dummy contrast data
agd_contrast <- agd_arm %>%
  group_by(studyn) %>%
  mutate(arm = 1:n(),
         ydiff = if_else(arm == 1, NA_real_, cont - first(cont)),
         ydiff_chr = if_else(arm == 1, NA_character_, trtc),
         ydiff_multi = if_else(arm %in% 1:2, NA_real_, cont - first(cont)),
         sediff = if_else(arm == 1, cont_pos, sqrt(cont_pos^2 + first(cont_pos)^2)),
         sediff_miss = if_else(arm == 1, NA_real_, sqrt(cont_pos^2 + first(cont_pos)^2))) %>%
  ungroup()

test_that("set_agd_contrast - continuous outcome checks work", {
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff), "Specify standard error")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", se = sediff), "Specify continuous outcome `y`")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff_chr, se = sediff), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = as.character(ydiff), se = sediff), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = trtc), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = as.character(sediff)), "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = cont_neg), "must be positive")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = -cont_pos), "must be positive")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = cont_inf), "cannot be infinite")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = cont_nan), "cannot be NaN")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = trtn, se = sediff), "without a specified baseline arm")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff_multi, se = sediff), "Multiple baseline arms")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = sediff_miss), "Standard error.+missing values on baseline arms")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = cbind(ydiff, ydiff), se = sediff), "must be a regular column")
  expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = cbind(sediff, sediff)), "must be a regular column")
  expect_equivalent(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff, se = sediff)$agd_contrast[, c(".y", ".se")],
                    transmute(agd_contrast, .y = ydiff, .se = sediff))
  expect_equivalent(set_agd_contrast(agd_contrast, "studyn", "trtc", y = "ydiff", se = "sediff")$agd_contrast[, c(".y", ".se")],
                    transmute(agd_contrast, .y = ydiff, .se = sediff))
  expect_equivalent(set_agd_contrast(agd_contrast, "studyn", "trtc", y = ydiff/2, se = sediff)$agd_contrast[, c(".y", ".se")],
                    transmute(agd_contrast, .y = ydiff/2, .se = sediff))
})

test_that("set_agd_contrast - sample size checks work", {
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = trtc),
               "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = as.character(trtn)),
               "must be numeric")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = cont),
               "must be integer")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = cont/2),
               "must be integer")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = disc_neg),
               "must be greater than zero")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = -disc),
               "must be greater than zero")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = disc_inf),
               "cannot be infinite")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = disc_nan),
               "cannot be NaN")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = cbind(disc, disc)),
               "must be a regular column")
  expect_error(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff, sample_size = list(disc)),
               "must be a regular column")
  expect_message(set_agd_contrast(agd_contrast, studyn, trtc, y = ydiff, se = sediff), "`sample_size` not provided")
})

test_that("set_agd_contrast - positive definite check", {
  agd_contrast_nonpd <- agd_contrast %>%
    group_by(studyn) %>%
    mutate(sediff = if_else(arm == 1, max(cont_pos)*2, sqrt(cont_pos^2 + first(cont_pos)^2))) %>%
    ungroup()

  expect_error(set_agd_contrast(agd_contrast_nonpd, studyn, trtc, y = ydiff, se = sediff),
               'not positive definite for study "2"')
  expect_error(set_agd_contrast(agd_contrast_nonpd, studyc, trtc, y = ydiff, se = sediff),
               'not positive definite for study "b"')

  agd_contrast_nonpd2 <- bind_rows(agd_contrast_nonpd,
                                   filter(agd_contrast_nonpd, studyn == 2) %>% mutate(studyn = 3, studyc = "c"))

  expect_error(set_agd_contrast(agd_contrast_nonpd2, studyn, trtc, y = ydiff, se = sediff),
               'not positive definite for studies "2" and "3"')
  expect_error(set_agd_contrast(agd_contrast_nonpd2, studyc, trtc, y = ydiff, se = sediff),
               'not positive definite for studies "b" and "c"')
})

test_that("set_* - take one and only one outcome", {
  m <- "specify one and only one outcome"
  expect_error(set_ipd(agd_arm, "studyn", "trtc", r = bin, y = cont), m)
  expect_error(set_agd_arm(agd_arm, "studyn", "trtc", r = disc, n = disc_p1, y = cont, se = cont_pos), m)
  # expect_error(set_agd_contrast(agd_contrast, "studyn", "trtc"), m)
})

# Reference treatment for survival outcomes tie-breaks by longest follow-up
bmax <- max(subset(agd_arm$cont_pos, agd_arm$trtf == "B")) > max(subset(agd_arm$cont_pos, agd_arm$trtf == "A"))

test_that("set_* `.trt` column is correct", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)

  expect_equal(set_ipd(agd_arm, studyc, 6, y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, 6, y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, 6, y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, 6, Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)

  expect_equal(set_ipd(agd_arm, studyc, "trtc", y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, "trtc", y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, "trtc", y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, "trtc", Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)

  expect_equal(set_ipd(agd_arm, studyc, factor(trtc), y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, factor(trtc), y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, factor(trtc), y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, factor(trtc), Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)


  expect_equal(set_ipd(agd_arm, studyc, trtf, y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf, y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyf, trtf, y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf, Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)

  expect_equal(set_ipd(agd_arm, studyc, trtf2, y = cont, trt_ref = "C")$ipd$.trt,
               agd_arm$trtf2)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf2, y = cont, se = cont_pos, trt_ref = "C")$agd_arm$.trt,
               agd_arm$trtf2)
  expect_equal(set_agd_contrast(agd_contrast, studyf, trtf2, y = ydiff, se = sediff, trt_ref = "C")$agd_contrast$.trt,
               agd_contrast$trtf2)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf2, Surv = Surv(cont_pos, bin), trt_ref = "C")$agd_arm$.trt,
               agd_arm$trtf2)

  # Check that unused factor levels are dropped
  expect_equal(set_ipd(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), y = cont)$ipd$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_arm(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), y = cont, se = cont_pos)$agd_arm$.trt,
               agd_arm$trtf)
  expect_equal(set_agd_contrast(agd_contrast, studyf, forcats::fct_expand(trtf, "zzz"), y = ydiff, se = sediff)$agd_contrast$.trt,
               agd_contrast$trtf)
  expect_equal(set_agd_surv(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), Surv = Surv(cont_pos, bin))$agd_arm$.trt,
               if (bmax) forcats::fct_relevel(agd_arm$trtf, "B", "A", "C") else agd_arm$trtf)

  expect_equal(set_ipd(agd_arm, studyc, forcats::fct_expand(trtf2, "zzz"), y = cont, trt_ref = "C")$ipd$.trt,
               agd_arm$trtf2)
  expect_equal(set_agd_arm(agd_arm, studyc, forcats::fct_expand(trtf2, "zzz"), y = cont, se = cont_pos, trt_ref = "C")$agd_arm$.trt,
               agd_arm$trtf2)
  expect_equal(set_agd_contrast(agd_contrast, studyf, forcats::fct_expand(trtf2, "zzz"), y = ydiff, se = sediff, trt_ref = "C")$agd_contrast$.trt,
               agd_contrast$trtf2)
  expect_equal(set_agd_surv(agd_arm, studyc, forcats::fct_expand(trtf2, "zzz"), Surv = Surv(cont_pos, bin), trt_ref = "C")$agd_arm$.trt,
               agd_arm$trtf2)
})

test_that("set_* `.study` column is correct", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, 2, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, 2, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, 2, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, 2, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, "studyc", trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, "studyc", trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, "studyc", trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, "studyc", trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, factor(studyc), trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, factor(studyc), trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, factor(studyc), trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, factor(studyc), trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, studyf, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, studyf, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, studyf, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, studyf, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, studyf2, trtc, y = cont)$ipd$.study,
               agd_arm$studyf2)
  expect_equal(set_agd_arm(agd_arm, studyf2, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf2)
  expect_equal(set_agd_contrast(agd_contrast, studyf2, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf2)
  expect_equal(set_agd_surv(agd_arm, studyf2, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf2)

  # Check that unused levels are dropped
  expect_equal(set_ipd(agd_arm, forcats::fct_expand(studyf, "zzz"), trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(agd_arm, forcats::fct_expand(studyf, "zzz"), trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(agd_contrast, forcats::fct_expand(studyf, "zzz"), trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(agd_arm, forcats::fct_expand(studyf, "zzz"), trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(agd_arm, forcats::fct_expand(studyf2, "zzz"), trtc, y = cont)$ipd$.study,
               agd_arm$studyf2)
  expect_equal(set_agd_arm(agd_arm, forcats::fct_expand(studyf2, "zzz"), trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf2)
  expect_equal(set_agd_contrast(agd_contrast, forcats::fct_expand(studyf2, "zzz"), trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf2)
  expect_equal(set_agd_surv(agd_arm, forcats::fct_expand(studyf2, "zzz"), trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf2)

  # Check reserved column names handled correctly
  aa <- mutate(agd_arm, .study = studyc)
  ac <- mutate(agd_contrast, .study = studyc)
  expect_equal(set_ipd(aa, .study, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(aa, .study, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(ac, .study, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(aa, .study, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_surv(aa, .study, trtc, Surv = Surv(cont_pos, bin), covariates = aa)$agd_arm$.study,
               agd_arm$studyf)

  expect_equal(set_ipd(aa, 26, trtc, y = cont)$ipd$.study,
               agd_arm$studyf)
  expect_equal(set_agd_arm(aa, 26, trtc, y = cont, se = cont_pos)$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_contrast(ac, 32, trtc, y = ydiff, se = sediff)$agd_contrast$.study,
               agd_contrast$studyf)
  expect_equal(set_agd_surv(aa, 26, trtc, Surv = Surv(cont_pos, bin))$agd_arm$.study,
               agd_arm$studyf)
  expect_equal(set_agd_surv(aa, 26, trtc, Surv = Surv(cont_pos, bin), covariates = aa)$agd_arm$.study,
               agd_arm$studyf)
})

test_that("set_* return default `treatments` factor", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$treatments,
               .default(factor(LETTERS[1:3])))
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$treatments,
               .default(factor(LETTERS[1:3])))
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff)$treatments,
               .default(factor(LETTERS[1:3])))
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin))$treatments,
               if (bmax) .default(factor(c("B", "A", "C"), levels = c("B", "A", "C"))) else .default(factor(LETTERS[1:3])))
})

test_that("set_* can set `trt_ref`", {
  f_BAC <- factor(LETTERS[c(2,1,3)], levels = LETTERS[c(2,1,3)])
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_ref = factor("B"))$treatments,
               f_BAC)

  # Using trtf sets original_levels attribute
  attr(f_BAC, "original_levels") <- LETTERS[1:3]
  expect_equal(set_ipd(agd_arm, studyc, trtf, y = cont, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf, y = cont, se = cont_pos, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtf, y = ydiff, se = sediff, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf, Surv = Surv(cont_pos, bin), trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_ipd(agd_arm, studyc, trtf, y = cont, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf, y = cont, se = cont_pos, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtf, y = ydiff, se = sediff, trt_ref = factor("B"))$treatments,
               f_BAC)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf, Surv = Surv(cont_pos, bin), trt_ref = factor("B"))$treatments,
               f_BAC)

  f_213 <- factor(c(2, 1, 3), levels = c(2, 1, 3))
  expect_equal(set_ipd(agd_arm, studyc, trtn, y = cont, trt_ref = 2)$treatments,
               f_213)
  expect_equal(set_agd_arm(agd_arm, studyc, trtn, y = cont, se = cont_pos, trt_ref = 2)$treatments,
               f_213)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtn, y = ydiff, se = sediff, trt_ref = 2)$treatments,
               f_213)
  expect_equal(set_agd_surv(agd_arm, studyc, trtn, Surv = Surv(cont_pos, bin), trt_ref = 2)$treatments,
               f_213)

  f_BCA <- factor(LETTERS[c(2,3,1)], levels = LETTERS[c(2,3,1)])
  attr(f_BCA, "original_levels") <- LETTERS[3:1]
  expect_equal(set_ipd(agd_arm, studyc, trtf2, y = cont, trt_ref = "B")$treatments,
               f_BCA)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf2, y = cont, se = cont_pos, trt_ref = "B")$treatments,
               f_BCA)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtf2, y = ydiff, se = sediff, trt_ref = "B")$treatments,
               f_BCA)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf2, Surv = Surv(cont_pos, bin), trt_ref = "B")$treatments,
               f_BCA)
  expect_equal(set_ipd(agd_arm, studyc, trtf2, y = cont, trt_ref = factor("B"))$treatments,
               f_BCA)
  expect_equal(set_agd_arm(agd_arm, studyc, trtf2, y = cont, se = cont_pos, trt_ref = factor("B"))$treatments,
               f_BCA)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtf2, y = ydiff, se = sediff, trt_ref = factor("B"))$treatments,
               f_BCA)
  expect_equal(set_agd_surv(agd_arm, studyc, trtf2, Surv = Surv(cont_pos, bin), trt_ref = factor("B"))$treatments,
               f_BCA)

  # Check that unused levels are dropped
  attr(f_BAC, "original_levels") <- c("A", "B", "C", "zzz")
  expect_equal(set_ipd(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), y = cont, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_arm(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), y = cont, se = cont_pos, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_contrast(agd_contrast, studyc, forcats::fct_expand(trtf, "zzz"), y = ydiff, se = sediff, trt_ref = "B")$treatments,
               f_BAC)
  expect_equal(set_agd_surv(agd_arm, studyc, forcats::fct_expand(trtf, "zzz"), Surv = Surv(cont_pos, bin), trt_ref = "B")$treatments,
               f_BAC)

  m <- "`trt_ref` does not match a treatment in the data.+Suitable values are:"
  expect_error(set_ipd(agd_arm, studyc, trtc, y = cont, trt_ref = 2), m)
  expect_error(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_ref = 2), m)
  expect_error(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_ref = 2), m)
  expect_error(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_ref = 2), m)
})

# Check classes when default reference treatment is not first in sort order
# Add new study to make B the default trt_ref
newstudy <- tibble(studyc = "c", trtc = c("B", "C"), tclassc = "b",
                   cont = rnorm(2), cont_pos = runif(2, 0, 1), bin = c(0, 1),
                   ydiff = c(NA, cont[2]), sediff = 1)
aa <- bind_rows(agd_arm, newstudy)
ac <- bind_rows(agd_contrast, newstudy)

test_that("set_* returns correct .trtclass column", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = tclassc)$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc)$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc)$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = tclassc)$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = tclassc))$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc))$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc))$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassc))$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)

  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = 11)$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = 11)$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = 11)$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = 11)$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = 11))$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = 11))$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = 11))$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = 11))$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)

  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = tclassf)$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassf)$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassf)$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = tclassf)$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = tclassf))$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = tclassf))$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = tclassf))$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassf))$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)

  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont, trt_ref = "B",
                       trt_class = tclassf2)$ipd$.trtclass,
               agd_arm$tclassf2)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_ref = "B",
                           trt_class = tclassf2)$agd_arm$.trtclass,
               agd_arm$tclassf2)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_ref = "B",
                                trt_class = tclassf2)$agd_contrast$.trtclass,
               agd_contrast$tclassf2)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_ref = "B",
                           trt_class = tclassf2)$agd_arm$.trtclass,
               agd_arm$tclassf2)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = tclassf2), trt_ref = "B")$ipd$.trtclass,
               agd_arm$tclassf2)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = tclassf2), trt_ref = "B")$agd_arm$.trtclass,
               agd_arm$tclassf2)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = tclassf2), trt_ref = "B")$agd_contrast$.trtclass,
               agd_contrast$tclassf2)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassf2), trt_ref = "B")$agd_arm$.trtclass,
               agd_arm$tclassf2)

  # Check that unused factor levels are dropped
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = forcats::fct_expand(tclassf, "zzz"))$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = forcats::fct_expand(tclassf, "zzz"))$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = forcats::fct_expand(tclassf, "zzz"))$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = forcats::fct_expand(tclassf, "zzz"))$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = forcats::fct_expand(tclassf, "zzz")))$ipd$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = forcats::fct_expand(tclassf, "zzz")))$agd_arm$.trtclass,
               agd_arm$tclassf)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = forcats::fct_expand(tclassf, "zzz")))$agd_contrast$.trtclass,
               agd_contrast$tclassf)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = forcats::fct_expand(tclassf, "zzz")))$agd_arm$.trtclass,
               if (bmax) forcats::fct_relevel(agd_arm$tclassf, "b") else agd_arm$tclassf)

  # Checks when default trt_ref not first in sort order
  expect_equal(set_ipd(aa, studyc, trtc, y = cont,
                       trt_class = tclassc)$ipd$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
  expect_equal(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc)$agd_arm$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
  expect_equal(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc)$agd_contrast$.trtclass,
               factor(ac$tclassc, levels = c("b", "a")))
  expect_equal(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = tclassc)$agd_arm$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
  expect_equal(combine_network(set_ipd(aa, studyc, trtc, y = cont,
                                       trt_class = tclassc))$ipd$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
  expect_equal(combine_network(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = tclassc))$agd_arm$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
  expect_equal(combine_network(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = tclassc))$agd_contrast$.trtclass,
               factor(ac$tclassc, levels = c("b", "a")))
  expect_equal(combine_network(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassc))$agd_arm$.trtclass,
               factor(aa$tclassc, levels = c("b", "a")))
})

test_that("set_* returns classes factor variable", {
  f_class <- factor(c("a", "b", "b"))
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont, trt_class = tclassc)$classes,
               f_class)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_class = tclassc)$classes,
               f_class)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_class = tclassc)$classes,
               f_class)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_class = tclassc)$classes,
               if (bmax) factor(c("b", "a", "b"), levels = c("b", "a")) else f_class)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont, trt_class = tclassc))$classes,
               f_class)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_class = tclassc))$classes,
               f_class)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff, trt_class = tclassc))$classes,
               f_class)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_class = tclassc))$classes,
               if (bmax) factor(c("b", "a", "b"), levels = c("b", "a")) else f_class)

  # Using tclassf sets original_levels attribute
  attr(f_class, "original_levels") <- c("a", "b")
  if (bmax) {
    f_classs <- factor(c("b", "a", "b"), levels = c("b", "a"))
    attr(f_classs, "original_levels") <- c("a", "b")
  }
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"))$classes,
               f_class)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"))$classes,
               f_class)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"))$classes,
               f_class)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"))$classes,
               if (bmax) f_classs else f_class)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")))$classes,
               if (bmax) f_classs else f_class)

  f_class2 <- factor(c("b", "a", "b"), levels = c("b", "a"))
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = tclassc, trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc, trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc, trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = tclassc, trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = tclassc), trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc), trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc), trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassc), trt_ref = "B")$classes,
               f_class2)

  # Checks when default trt_ref not first in sort order
  expect_equal(set_ipd(aa, studyc, trtc, y = cont,
                       trt_class = tclassc)$classes,
               f_class2)
  expect_equal(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = tclassc)$classes,
               f_class2)
  expect_equal(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = tclassc)$classes,
               f_class2)
  expect_equal(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = tclassc)$classes,
               f_class2)
  expect_equal(combine_network(set_ipd(aa, studyc, trtc, y = cont,
                                       trt_class = tclassc))$classes,
               f_class2)
  expect_equal(combine_network(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = tclassc))$classes,
               f_class2)
  expect_equal(combine_network(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = tclassc))$classes,
               f_class2)
  expect_equal(combine_network(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = tclassc))$classes,
               f_class2)

  attr(f_class2, "original_levels") <- c("a", "b")
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                       trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                           trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                                trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                           trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)

  # Checks when default trt_ref not first in sort order
  expect_equal(set_ipd(aa, studyc, trtc, y = cont,
                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                       trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                           trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                                trt_ref = "B")$classes,
               f_class2)
  expect_equal(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b"),
                           trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_ipd(aa, studyc, trtc, y = cont,
                                       trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_arm(aa, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_contrast(ac, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)
  expect_equal(combine_network(set_agd_surv(aa, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = recode_factor(trtc, A = "a", B = "b", C = "b")),
                               trt_ref = "B")$classes,
               f_class2)

  f_class3 <- factor(c("B", "A", "C"), levels = c("B", "C", "A"))
  attr(f_class3, "original_levels") <- c("C", "B", "A")
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = trtf2,
                       trt_ref = "B")$classes,
               f_class3)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = trtf2,
                           trt_ref = "B")$classes,
               f_class3)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = trtf2,
                                trt_ref = "B")$classes,
               f_class3)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = trtf2,
                           trt_ref = "B")$classes,
               f_class3)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = trtf2),
                               trt_ref = "B")$classes,
               f_class3)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = trtf2),
                               trt_ref = "B")$classes,
               f_class3)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = trtf2),
                               trt_ref = "B")$classes,
               f_class3)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = trtf2),
                               trt_ref = "B")$classes,
               f_class3)

  # Check that unused levels are dropped
  attr(f_class, "original_levels") <- c("a", "b", "zzz")
  if (bmax) attr(f_classs, "original_levels") <- c("a", "b", "zzz")
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont,
                       trt_class = forcats::fct_expand(tclassf, "zzz"))$classes,
               f_class)
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                           trt_class = forcats::fct_expand(tclassf, "zzz"))$classes,
               f_class)
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                trt_class = forcats::fct_expand(tclassf, "zzz"))$classes,
               f_class)
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                           trt_class = forcats::fct_expand(tclassf, "zzz"))$classes,
               if (bmax) f_classs else f_class)
  expect_equal(combine_network(set_ipd(agd_arm, studyc, trtc, y = cont,
                                       trt_class = forcats::fct_expand(tclassf, "zzz")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos,
                                           trt_class = forcats::fct_expand(tclassf, "zzz")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff,
                                                trt_class = forcats::fct_expand(tclassf, "zzz")))$classes,
               f_class)
  expect_equal(combine_network(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin),
                                           trt_class = forcats::fct_expand(tclassf, "zzz")))$classes,
               if (bmax) f_classs else f_class)
})

test_that("set_* checks for bad class variable work", {
  aa2 <- agd_arm
  aa2$tclassn <- c(1, 2, 2, 1, 2) # Trt 2 and 3 in two classes
  aa2$tclassc[1] <- NA

  ac2 <- agd_contrast
  ac2$tclassn <- c(1, 2, 2, 1, 2) # Trt 2 and 3 in two classes
  ac2$tclassc[1] <- NA

  m <- "Treatment present in more than one class"

  expect_error(set_ipd(aa2, studyc, trtc, y = cont, trt_class = tclassn), m)
  expect_error(set_agd_arm(aa2, studyc, trtc, y = cont, se = cont_pos, trt_class = tclassn), m)
  expect_error(set_agd_contrast(ac2, studyc, trtc, y = ydiff, se = sediff, trt_class = tclassn), m)
  expect_error(set_agd_surv(aa2, studyc, trtc, Surv = Surv(cont_pos, bin), trt_class = tclassn), m)

  m2 <- "cannot contain missing values"

  expect_error(set_ipd(aa2, studyc, trtc, y = cont, trt_class = tclassc), m2)
  expect_error(set_agd_arm(aa2, studyc, trtc, y = cont, se = cont_pos, trt_class = tclassc), m2)
  expect_error(set_agd_contrast(ac2, studyc, trtc, y = ydiff, se = sediff, trt_class = tclassc), m2)
  expect_error(set_agd_surv(aa2, studyc, trtc, Surv = Surv(cont_pos, bin), trt_class = tclassc), m2)

  m3 <- "must be a regular column"
  expect_error(set_ipd(agd_arm, studyc, trtc, y = cont, trt_class = cbind(trtc, trtc)), m3)
  expect_error(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos, trt_class = cbind(trtc, trtc)), m3)
  expect_error(set_agd_contrast(agd_contrast, studyc, trtc, y = cont, se = cont_pos, trt_class = cbind(trtc, trtc)), m3)
  expect_error(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin), trt_class = cbind(trtc, trtc)), m3)
})

test_that("set_* return `studies` factor", {
  expect_equal(set_ipd(agd_arm, studyc, trtc, y = cont)$studies,
               factor(letters[1:2]))
  expect_equal(set_agd_arm(agd_arm, studyc, trtc, y = cont, se = cont_pos)$studies,
               factor(letters[1:2]))
  expect_equal(set_agd_contrast(agd_contrast, studyc, trtc, y = ydiff, se = sediff)$studies,
               factor(letters[1:2]))
  expect_equal(set_agd_surv(agd_arm, studyc, trtc, Surv = Surv(cont_pos, bin))$studies,
               factor(letters[1:2]))
})

make_na <- function(x, n) {
  x[sample.int(length(x), n)] <- NA
  return(x)
}

test_that("set_* error if outcomes contain missing values", {
  agd_arm_miss <- agd_arm %>% mutate_at(vars(cont:bin), ~make_na(., 1))

  m <- "contains missing values"

  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, y = cont, se = cont_pos), m)
  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, r = disc, n = disc_p1), m)
  expect_error(set_agd_arm(agd_arm_miss, studyn, trtn, r = disc, E = cont_pos), m)

  expect_error(set_ipd(agd_arm_miss, studyn, trtn, y = cont), m)
  expect_error(set_ipd(agd_arm_miss, studyn, trtn, r = disc), m)
  expect_error(set_ipd(agd_arm_miss, studyn, trtn, r = disc, E = cont_pos), m)

  # For set_agd_contrast, multiple missing y will trigger multiple baselines error
  # Missing se on baseline arm will trigger covariance error if >2 arms (otherwise fine)
  agd_contrast_miss <- agd_contrast %>% mutate(sediff = if_else(arm > 1, NA_real_, sediff))
  expect_error(set_agd_contrast(agd_contrast_miss, studyn, trtn, y = ydiff, se = sediff), m)
})

