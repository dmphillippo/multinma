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

test_that("treatments must be in network", {
  m1 <- "`trt1` does not match a treatment in the network"
  m2 <- "`trt2` does not match a treatment in the network"

  expect_error(has_direct(thrombo_net, "bad", "SK"), m1)
  expect_error(has_direct(thrombo_net, "SK", "bad"), m2)
  expect_error(has_indirect(thrombo_net, "bad", "SK"), m1)
  expect_error(has_indirect(thrombo_net, "SK", "bad"), m2)
})

test_that("trt1 and trt2 must be different", {
  m <- "`trt1` and `trt2` cannot be the same treatment"

  expect_error(has_direct(thrombo_net, "SK", "SK"), m)
  expect_error(has_indirect(thrombo_net, "SK", "SK"), m)
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

test_that("get_nodesplits() produces correct output", {
  # Output from gemtc::mtc.nodesplit.comparisons() on thrombo network
  ns_gemtc <- tibble::tribble(
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

  expect_identical(get_nodesplits(thrombo_net2), ns_gemtc)
})
