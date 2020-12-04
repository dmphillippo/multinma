library(multinma)
library(dplyr)

test_that("combine_network error if not passed nma_data objects", {
  msg <- "Expecting to combine objects of class `nma_data`, created using set_\\* functions"

  expect_error(combine_network(1), msg)
  expect_error(combine_network(1, 2), msg)
  expect_error(combine_network(1, set_agd_arm(smoking, studyn, trtc, r = r, n = n)), msg)
  expect_error(combine_network(set_agd_arm(smoking, studyn, trtc, r = r, n = n), 1), msg)
})

# Dummy data
agd_arm <- tibble(
  studyn = c(1, 1, 2, 2, 2),
  studyc = letters[studyn],
  studyf = factor(studyc),
  trtn = c(1, 2, 1, 2, 3),
  trtc = LETTERS[trtn],
  trtf = factor(trtc),
  tclassn = c(1, 2, 1, 2, 2),
  tclassc = letters[tclassn],
  tclassf = factor(tclassc, levels = letters[1:3]),
  y = rnorm(5),
  se = runif(5)
)
net_a_a <- set_agd_arm(agd_arm, studyf, trtf, y = y, se = se)

agd_contrast <- tibble(
  studyn = c(3, 3, 3),
  studyc = letters[studyn],
  studyf = factor(studyc),
  trtn = c(2, 3, 4),
  trtc = LETTERS[trtn],
  trtf = factor(trtc),
  tclassn = c(2, 2, 3),
  tclassc = letters[tclassn],
  tclassf = factor(tclassc, levels = letters[1:3]),
  tclassn_bad = c(3, 3, 3),
  y = c(NA, rnorm(2)),
  se = runif(3)
)
net_a_c <- set_agd_contrast(agd_contrast, studyf, trtf, y = y, se = se)

ipd <- tibble(
  studyn = c(4, 4, 4, 5, 5),
  studyc = letters[studyn],
  studyf = factor(studyc),
  trtn = c(1, 2, 3, 1, 4),
  trtc = LETTERS[trtn],
  trtf = factor(trtc),
  tclassn = c(1, 2, 2, 1, 3),
  tclassc = letters[tclassn],
  tclassf = factor(tclassc, levels = letters[1:3]),
  tclassn_bad = c(1, 3, 3, 3, 4),
  y = rnorm(5)
)
net_i <- set_ipd(ipd, studyf, trtf, y = y)

test_that("combine_network produces combined treatment, class, and study factors", {
  c1 <- combine_network(net_a_a, net_i)
  expect_equal(c1$treatments, .default(factor(LETTERS[1:4])))
  expect_equal(levels(c1$agd_arm$.trt), LETTERS[1:4])
  expect_equal(levels(c1$ipd$.trt), LETTERS[1:4])
  expect_equal(c1$studies, factor(letters[c(1, 2, 4, 5)]))
  expect_equal(levels(c1$agd_arm$.study), letters[c(1, 2, 4, 5)])
  expect_equal(levels(c1$ipd$.study), letters[c(1, 2, 4, 5)])

  c2 <- combine_network(net_a_a, net_i, net_a_c)
  expect_equal(c2$treatments, .default(factor(LETTERS[c(2, 1, 3, 4)],
                                              levels = LETTERS[c(2, 1, 3, 4)])))
  expect_equal(levels(c2$agd_arm$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(levels(c2$agd_contrast$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(levels(c2$ipd$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(c2$studies, factor(letters[1:5]))
  expect_equal(levels(c2$agd_arm$.study), letters[1:5])
  expect_equal(levels(c2$agd_contrast$.study), letters[1:5])
  expect_equal(levels(c2$ipd$.study), letters[1:5])

  c1_classed <- combine_network(
    set_agd_arm(agd_arm, studyf, trtf, y = y, se = se, trt_class = tclassc),
    set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassc)
  )
  # Reference trt is A
  expect_equal(c1_classed$treatments, .default(factor(LETTERS[1:4])))
  expect_equal(levels(c1_classed$agd_arm$.trt), LETTERS[1:4])
  expect_equal(levels(c1_classed$ipd$.trt), LETTERS[1:4])
  expect_equal(c1_classed$studies, factor(letters[c(1, 2, 4, 5)]))
  expect_equal(levels(c1_classed$agd_arm$.study), letters[c(1, 2, 4, 5)])
  expect_equal(levels(c1_classed$ipd$.study), letters[c(1, 2, 4, 5)])
  expect_equal(c1_classed$classes, factor(letters[c(1, 2, 2, 3)]))
  expect_equal(levels(c1_classed$agd_arm$.trtclass), letters[1:3])
  expect_equal(levels(c1_classed$ipd$.trtclass), letters[1:3])

  c2_classed <- combine_network(
    set_agd_arm(agd_arm, studyf, trtf, y = y, se = se, trt_class = tclassc),
    set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassc),
    set_agd_contrast(agd_contrast, studyf, trtf, y = y, se = se, trt_class = tclassc)
  )
  # Reference treatment is B
  expect_equal(c2_classed$treatments, .default(factor(LETTERS[c(2, 1, 3, 4)],
                                                      levels = LETTERS[c(2, 1, 3, 4)])))
  expect_equal(levels(c2_classed$agd_arm$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(levels(c2_classed$agd_contrast$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(levels(c2_classed$ipd$.trt), LETTERS[c(2, 1, 3, 4)])
  expect_equal(c2_classed$studies, factor(letters[1:5]))
  expect_equal(levels(c2_classed$agd_arm$.study), letters[1:5])
  expect_equal(levels(c2_classed$agd_contrast$.study), letters[1:5])
  expect_equal(levels(c2_classed$ipd$.study), letters[1:5])
  expect_equal(c2_classed$classes, factor(letters[c(2, 1, 2, 3)], levels = c(letters[c(2, 1, 3)])))
  expect_equal(levels(c2_classed$agd_arm$.trtclass), letters[c(2, 1, 3)])
  expect_equal(levels(c2_classed$ipd$.trtclass), letters[c(2, 1, 3)])
  expect_equal(levels(c2_classed$agd_contrast$.trtclass), letters[c(2, 1, 3)])
})

test_that("combine_network can set alternative trt_ref", {
  c1 <- combine_network(net_a_a, net_i, net_a_c, trt_ref = "C")
  expect_equal(c1$treatments, factor(LETTERS[c(3, 1, 2, 4)], levels = LETTERS[c(3, 1, 2, 4)]))
  expect_equal(levels(c1$agd_arm$.trt), LETTERS[c(3, 1, 2, 4)])
  expect_equal(levels(c1$agd_contrast$.trt), LETTERS[c(3, 1, 2, 4)])
  expect_equal(levels(c1$ipd$.trt), LETTERS[c(3, 1, 2, 4)])

  expect_error(combine_network(net_a_a, net_i, net_a_c, trt_ref = 2),
               "does not match a treatment.*Suitable values are: A, B, C, D")
  expect_error(combine_network(net_a_a, net_a_c,
                               set_ipd(mutate(ipd, trtf = factor(LETTERS[3:7])),
                                       studyf, trtf, y = y),
                               trt_ref = 2),
               "does not match a treatment.*Suitable values are: A, B, C, D, E, \\.\\.\\.")

  expect_error(combine_network(net_a_a, net_a_c,
                               set_ipd(mutate(ipd, studyf = factor("a")),
                                       studyf, trtf, y = y)),
               "Studies with same label found in multiple data sources: a")
  expect_error(combine_network(net_a_a, net_a_a),
               "Studies with same label found in multiple data sources: a, b")

  c1_classed <- combine_network(
    set_agd_arm(agd_arm, studyf, trtf, y = y, se = se, trt_class = tclassc),
    set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassc),
    set_agd_contrast(agd_contrast, studyf, trtf, y = y, se = se, trt_class = tclassc),
    trt_ref = "C"
  )
  expect_equal(c1_classed$treatments, factor(LETTERS[c(3, 1, 2, 4)], levels = LETTERS[c(3, 1, 2, 4)]))
  expect_equal(levels(c1_classed$agd_arm$.trt), LETTERS[c(3, 1, 2, 4)])
  expect_equal(levels(c1_classed$agd_contrast$.trt), LETTERS[c(3, 1, 2, 4)])
  expect_equal(levels(c1_classed$ipd$.trt), LETTERS[c(3, 1, 2, 4)])
  expect_equal(c1_classed$classes, factor(letters[c(2, 1, 2, 3)], levels = letters[c(2, 1, 3)]))
  expect_equal(levels(c1_classed$agd_arm$.trtclass), letters[c(2, 1, 3)])
  expect_equal(levels(c1_classed$agd_contrast$.trtclass), letters[c(2, 1, 3)])
  expect_equal(levels(c1_classed$ipd$.trtclass), letters[c(2, 1, 3)])
})

test_that("combine_network error if outcomes do not match for same data source type", {
  m <- "Multiple outcome types present"
  dat_a <- tibble(study = 1, trt = 1:2, r = 1, n = 1, y = 1, se = 1)
  dat_b <- tibble(study = 2, trt = 2:3, r = 1, n = 1, y = 1, se = 1)

  expect_error(combine_network(set_ipd(dat_a, study, trt, r = r),
                               set_ipd(dat_b, study, trt, y = y)), m)
  expect_error(combine_network(set_agd_arm(dat_a, study, trt, r = r, n = n),
                               set_agd_arm(dat_b, study, trt, y = y, se = se)), m)
})

test_that("combine_network error if mismatch outcomes across data types", {
  m <- "Combining.+not supported"
  dat_a <- tibble(study = 1, trt = 1:2, r = 1, n = 1, E = 1, y = 1, se = 1)
  dat_b <- tibble(study = 2, trt = 2:3, r = 1, n = 1, E = 1, y = 1, se = 1)

  expect_error(combine_network(set_ipd(dat_a, study, trt, r = r),
                               set_agd_arm(dat_b, study, trt, y = y, se = se)), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(n, r, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, y = y, se = se)), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(n, r, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = r, n = n)), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = r, E = E),
                               set_agd_arm(dat_b, study, trt, y = y, se = se)), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, y = y),
                               set_agd_arm(dat_b, study, trt, r = r, n = n)), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, y = y),
                               set_agd_arm(dat_b, study, trt, r = multi(n, r, inclusive = TRUE))), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = r),
                               set_agd_arm(dat_b, study, trt, r = multi(n, r, inclusive = TRUE))), m)
  expect_error(combine_network(set_ipd(dat_a, study, trt, y = y),
                               set_agd_arm(dat_b, study, trt, r = r, E = E)), m)
})

test_that("combine_network error if multinomial outcomes mismatched", {
  dat_a <- tibble(study = 1, trt = 1:2, a = 1, b = 1, c = 0)
  dat_b <- tibble(study = 2, trt = 2:3, a = 1, b = 1, c = 0)

  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(a, b, inclusive = TRUE))),
               "different numbers of categories")
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(b, a, c, inclusive = TRUE))),
               "different category labels")
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(A = a, B = b, C = c, inclusive = TRUE))),
               "different category labels")

  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_ipd(dat_b, study, trt, r = multi(a, b, inclusive = TRUE))),
               "different numbers of categories")
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_ipd(dat_b, study, trt, r = multi(b, a, c, inclusive = TRUE))),
               "different category labels")
  expect_error(combine_network(set_ipd(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_ipd(dat_b, study, trt, r = multi(A = a, B = b, C = c, inclusive = TRUE))),
               "different category labels")

  expect_error(combine_network(set_agd_arm(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(a, b, inclusive = TRUE))),
               "different numbers of categories")
  expect_error(combine_network(set_agd_arm(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(b, a, c, inclusive = TRUE))),
               "different category labels")
  expect_error(combine_network(set_agd_arm(dat_a, study, trt, r = multi(a, b, c, inclusive = TRUE)),
                               set_agd_arm(dat_b, study, trt, r = multi(A = a, B = b, C = c, inclusive = TRUE))),
               "different category labels")

  # Check combining ordered and competing outcomes is disallowed, when competing are implemented
})

test_that("combine_network error if treatment classes do not match across sources", {
  m <- "Treatment present in more than one class"

  expect_error(combine_network(
    set_agd_arm(agd_arm, studyf, trtf, y = y, se = se, trt_class = tclassn),
    set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassn_bad)), m)

  expect_error(combine_network(
    set_agd_arm(agd_arm, studyf, trtf, y = y, se = se, trt_class = tclassn),
    set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassn_bad),
    set_agd_contrast(agd_contrast, studyf, trtf, y = y, se = se, trt_class = tclassn_bad)), m)
})

test_that("combine_network warns if not all sources have treatment classes", {
  expect_warning(combine_network(net_a_a,
                                 set_ipd(ipd, studyf, trtf, y = y, trt_class = tclassc)),
                 "Not all data sources have defined treatment classes")
})
