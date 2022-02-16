library(multinma)
library(dplyr)

test_that("network must be a network", {
  m <- "`network` must be an `nma_data` object"
  expect_error(has_direct(list()), m)
  expect_error(has_indirect(list()), m)
  expect_error(get_nodesplits(list()), m)
})

thrombo_net <- set_agd_arm(thrombolytics, studyn, trtc,
                           r = r, n = n)
thrombo_net2 <- set_agd_arm(thrombolytics, studyn, trtn,
                            r = r, n = n)

test_that("treatments must be scalars", {
  m <- "should be a single integer, string, or factor, naming a treatment"

  expect_error(has_direct(thrombo_net, 1:2, "SK"), m)
  expect_error(has_direct(thrombo_net, "SK", 1:2), m)
  expect_error(has_indirect(thrombo_net, 1:2, "SK"), m)
  expect_error(has_indirect(thrombo_net, "SK", 1:2), m)
})

test_that("treatments must be in network", {
  m1 <- "`trt1` does not match a treatment in the network"
  m2 <- "`trt2` does not match a treatment in the network"

  expect_error(has_direct(thrombo_net, "bad", "SK"), m1)
  expect_error(has_direct(thrombo_net, "SK", "bad"), m2)
  expect_error(has_indirect(thrombo_net, "bad", "SK"), m1)
  expect_error(has_indirect(thrombo_net, "SK", "bad"), m2)

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = c("bad", "SK"),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    "The `nodesplit` treatment comparison should match two treatments in the network")
})

test_that("trt1 and trt2 must be different", {
  m <- "`trt1` and `trt2` cannot be the same treatment"

  expect_error(has_direct(thrombo_net, "SK", "SK"), m)
  expect_error(has_indirect(thrombo_net, "SK", "SK"), m)

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = c("SK", "SK"),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    "`nodesplit` comparison cannot be the same treatment against itself")
})

test_that("has_direct and has_indirect outputs are correct", {
  expect_identical(has_direct(thrombo_net, "TNK", "Acc t-PA"), TRUE)
  expect_identical(has_direct(thrombo_net, "TNK", "t-PA"), FALSE)
  expect_identical(has_indirect(thrombo_net, "TNK", "Acc t-PA"), FALSE)
  expect_identical(has_indirect(thrombo_net, "TNK", "t-PA"), TRUE)

  expect_identical(has_direct(thrombo_net2, 6, 3), TRUE)
  expect_identical(has_direct(thrombo_net2, 6, 2), FALSE)
  expect_identical(has_indirect(thrombo_net2, 6, 3), FALSE)
  expect_identical(has_indirect(thrombo_net2, 6, 2), TRUE)
})

# Output from gemtc::mtc.nodesplit.comparisons() on thrombo network
ns_thrombo_gemtc <- tibble::tribble(
  ~trt1, ~trt2,
  1,     2,
  1,     3,
  1,     5,
  1,     7,
  1,     8,
  1,     9,
  2,     7,
  2,     8,
  2,     9,
  3,     4,
  3,     5,
  3,     7,
  3,     8,
  3,     9
) %>%
  mutate(trt1 = factor(trt1, levels = levels(thrombo_net2$treatments)),
         trt2 = factor(trt2, levels = levels(thrombo_net2$treatments)))

test_that("get_nodesplits() produces correct output for thombolytics network", {
  expect_identical(get_nodesplits(thrombo_net2), ns_thrombo_gemtc)
})

park_net <- set_agd_arm(parkinsons, studyn, trtn, y = y, se = se, trt_ref = 1)

# Compare to results in van Valkenhoef paper
ns_park_vv <- tibble::tribble(
  ~trt1, ~trt2,
  1, 3,
  1, 4,
  2, 4,
  3, 4
) %>%
  mutate(trt1 = factor(trt1, levels = levels(park_net$treatments)),
         trt2 = factor(trt2, levels = levels(park_net$treatments)))

test_that("get_nodesplits() produces correct output for parkinsons network", {
  expect_identical(get_nodesplits(park_net), ns_park_vv)
})

test_that("get_nodesplits() handles repeated treatment arms", {
  thrombo_net_rep <- set_agd_arm(rbind(thrombolytics, thrombolytics[c(2, 5), ]),
                                 studyn, trtn,
                                 r = r, n = n)

  expect_identical(get_nodesplits(thrombo_net_rep), ns_thrombo_gemtc)
})

onestudy <- data.frame(study = 1, trt = 1:3, r = 1, n = 1)
pair_net <- set_agd_arm(onestudy[1:2, ], study, trt, r = r, n = n)
multi_net <- set_agd_arm(onestudy, study, trt, r = r, n = n)

test_that("get_nodesplits() returns an empty tibble if no splits to be done", {
  expect_identical(get_nodesplits(pair_net),
                   tibble(trt1 = factor(levels = levels(pair_net$treatments)),
                          trt2 = factor(levels = levels(pair_net$treatments))))

  expect_identical(get_nodesplits(multi_net),
                   tibble(trt1 = factor(levels = levels(multi_net$treatments)),
                          trt2 = factor(levels = levels(multi_net$treatments))))
})

test_that("has_direct and has_indirect work with one study / pairwise MA", {
  expect_identical(has_direct(pair_net, 1, 2), TRUE)
  expect_identical(has_direct(multi_net, 1, 2), TRUE)

  expect_identical(has_indirect(pair_net, 1, 2), FALSE)
  expect_identical(has_indirect(multi_net, 1, 2), FALSE)
})

test_that("nma() nodesplit error if can't split given comparison", {
  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = c("TNK", "t-PA"),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    "no direct evidence"
  )

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = c("TNK", "Acc t-PA"),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    "no independent indirect evidence"
  )
})

test_that("nma() nodesplit error if no comparisons to split", {
  m <- "No comparisons to node-split"

  expect_error(
    nma(pair_net,
        consistency = "nodesplit",
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m
  )

  expect_error(
    nma(multi_net,
        consistency = "nodesplit",
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m
  )
})

test_that("nma() nodesplit warning about ignoring comparisons", {
  expect_warning(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = data.frame(trt1 = "Acc t-PA", trt2 = c("TNK", "ASPAC")),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10),
        test_grad = TRUE),
    "Ignoring node-split comparisons.+Acc t-PA vs\\. TNK"
  )
})

test_that("nma() nodesplit argument validation", {
  m1 <- "The data frame passed to `nodesplit` should have two columns"

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = data.frame(),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m1
  )

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = data.frame(a = 1, b = 2, c = 3),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m1
  )

  m2 <- "`nodesplit` should either be a length 2 vector or a 2 column data frame"

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = list(),
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m2
  )

  expect_error(
    nma(thrombo_net,
        consistency = "nodesplit",
        nodesplit = 1:3,
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10)),
    m2
  )
})

thrombo_ns1 <- nma(thrombo_net,
                   consistency = "nodesplit",
                   nodesplit = c("SK", "t-PA"),
                   prior_intercept = normal(scale = 10),
                   prior_trt = normal(scale = 10),
                   test_grad = TRUE)

thrombo_ns2<- nma(thrombo_net,
                  consistency = "nodesplit",
                  nodesplit = data.frame(trt1 = "SK", trt2 = "t-PA"),
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  test_grad = TRUE)

test_that("nma() nodesplit produces correct classes", {
  expect_s3_class(thrombo_ns1, "nma_nodesplit")
  expect_s3_class(thrombo_ns2, "nma_nodesplit_df")
})

test_that("summary method consistency argument requires consistency model", {
  m <- "should be a fitted consistency model"
  expect_error(summary(thrombo_ns1, consistency = thrombo_ns1), m)
  expect_error(summary(thrombo_ns2, consistency = thrombo_ns1), m)
  expect_error(summary(thrombo_ns1, consistency = "bad"), m)
  expect_error(summary(thrombo_ns2, consistency = "bad"), m)
})

test_that("summary method checks nodesplit and consistency models are compatible", {
  skip_on_cran()

  m <- "does not match the node-splitting model"

  suppressWarnings({
    thrombo_ns1 <- nma(thrombo_net,
                       consistency = "nodesplit",
                       nodesplit = c("SK", "t-PA"),
                       prior_intercept = normal(scale = 10),
                       prior_trt = normal(scale = 10),
                       iter = 1)

    thrombo_ns2<- nma(thrombo_net,
                      consistency = "nodesplit",
                      nodesplit = data.frame(trt1 = "SK", trt2 = "t-PA"),
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      iter = 1)

    thrombo_re <- nma(thrombo_net,
                      trt_effects = "random",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_het = half_normal(2),
                      iter = 1)

    thrombo_reg <- nma(thrombo_net,
                       regression = ~I(rnorm(102)),
                       prior_intercept = normal(scale = 10),
                       prior_trt = normal(scale = 10),
                       prior_reg = normal(scale = 10),
                       iter = 1)
  })

  expect_error(summary(thrombo_ns1, consistency = thrombo_re), m)
  expect_error(summary(thrombo_ns2, consistency = thrombo_re), m)
  expect_error(summary(thrombo_ns1, consistency = thrombo_reg), m)
  expect_error(summary(thrombo_ns2, consistency = thrombo_reg), m)
})
