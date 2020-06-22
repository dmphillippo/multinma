
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

test_that("Error if not ML-NMR model", {
  expect_error(plot_integration_error(smk_fit_RE),
               "Expecting a `stan_mlnmr` object, created by fitting a ML-NMR model with numerical integration")
})
