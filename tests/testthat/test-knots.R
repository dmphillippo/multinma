
ndmm_net <- combine_network(
  set_ipd(ndmm_ipd, study = study, trt = trt,
          Surv = Surv(eventtime, status)),
  set_agd_surv(ndmm_agd, study = study, trt = trt,
               Surv = Surv(eventtime, status)),
  set_agd_contrast(data.frame(study = "A", trt = c("Pbo", "Len"), y = c(NA, 0), se = c(NA, 1), ss = 10),
                   study = study, trt = trt, y = y, se = se, sample_size = ss)
)

test_that("make_knots()", {

  expect_error(make_knots("test"), "must be an `nma_data` object")
  expect_error(make_knots(ndmm_net, n_knots = "1"), "`n_knots` must be a positive integer")
  expect_error(make_knots(ndmm_net, type = "wrong"), "`type` must be one of")

  # defaults
  k <- make_knots(ndmm_net)
  expect_named(k, setdiff(levels(ndmm_net$studies), "A"))

  attal <- dplyr::filter(ndmm_ipd, study == "Attal2012")
  attal_obs <- dplyr::filter(attal, status == 1)

  expect_identical(k$Attal2012,
                   c(0,
                     quantile(attal_obs$eventtime, probs = (1:7)/8, names = FALSE),
                     max(attal$eventtime)))

  # setting n_knots
  k2 <- make_knots(ndmm_net, n_knots = 3)

  expect_named(k2, setdiff(levels(ndmm_net$studies), "A"))
  expect_identical(k2$Attal2012,
                   c(0,
                     quantile(attal_obs$eventtime, probs = (1:3)/4, names = FALSE),
                     max(attal$eventtime)))

  # common knots
  k3 <- make_knots(ndmm_net, type = "quantile_common")
  expect_named(k3, setdiff(levels(ndmm_net$studies), "A"))
  expect_all_equal(k3, k3[1])
})

test_that("knots()", {
  ndmm_fit_mspline <- suppressWarnings(nma(ndmm_net,
                                           likelihood = "mspline",
                                           prior_intercept = normal(0, 100),
                                           prior_trt = normal(0, 10),
                                           prior_aux = half_normal(1),
                                           test_grad = TRUE))

  expect_identical(knots(ndmm_fit_mspline), make_knots(ndmm_net))
  expect_identical(knots(ndmm_fit_mspline, type = "internal"),
                   purrr::map(make_knots(ndmm_net), ~.[-c(1, 9)]))
  expect_identical(knots(ndmm_fit_mspline, type = "boundary"),
                   purrr::map(make_knots(ndmm_net), ~.[c(1, 9)]))

  ndmm_fit_mspline2 <- suppressWarnings(nma(ndmm_net,
                                           likelihood = "mspline",
                                           n_knots = 3,
                                           prior_intercept = normal(0, 100),
                                           prior_trt = normal(0, 10),
                                           prior_aux = half_normal(1),
                                           test_grad = TRUE))

  expect_identical(knots(ndmm_fit_mspline2), make_knots(ndmm_net, n_knots = 3))

  ndmm_fit_mspline_nphr <- suppressWarnings(nma(ndmm_net,
                                                likelihood = "mspline",
                                                aux_regression = ~.trt,
                                                prior_intercept = normal(0, 100),
                                                prior_trt = normal(0, 10),
                                                prior_aux = half_normal(1),
                                                prior_aux_reg = normal(0, 10),
                                                test_grad = TRUE))

  expect_identical(knots(ndmm_fit_mspline_nphr), make_knots(ndmm_net, type = "quantile_common"))

  ndmm_fit_pexp <- suppressWarnings(nma(ndmm_net,
                                           likelihood = "pexp",
                                           prior_intercept = normal(0, 100),
                                           prior_trt = normal(0, 10),
                                           prior_aux = half_normal(1),
                                           test_grad = TRUE))

  expect_identical(knots(ndmm_fit_pexp), make_knots(ndmm_net))

  ndmm_fit_weib <- suppressWarnings(nma(ndmm_net,
                                           likelihood = "weibull",
                                           prior_intercept = normal(0, 100),
                                           prior_trt = normal(0, 10),
                                           prior_aux = half_normal(1),
                                           test_grad = TRUE))

  expect_error(knots(ndmm_fit_weib), "No knots present")
})

