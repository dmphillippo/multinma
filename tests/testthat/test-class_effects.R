sa_net <- set_agd_contrast(social_anxiety, studyc,
                           trtc,
                           y = y,
                           se = se,
                           trt_class = classc,
                           trt_ref = "Waitlist")

sa_fit_EXclass_RE <- suppressWarnings(nma(sa_net, trt_effects = "random",
                                          prior_trt = normal(0, 100),
                                          prior_het = half_normal(5),
                                          class_effects = "exchangeable",
                                          iter = 10,
                                          prior_class_sd = normal(0.33,0.1),
                                          class_sd =
                                            list(`Exercise and SH no support` =
                                                   c("Exercise promotion", "Self-help no support"),
                                                 `SSRIs and NSSA` =
                                                   c("SSRI/SNRI", "NSSA"),
                                                 `Psychodynamic & Other psychological therapies` =
                                                   c("Psychodynamic psychotherapy", "Other psychological therapies"))))

test_that("class_mean and class_sd parameters are named correctly", {
  # Extract summaries for class_mean
  cm_summary <- tibble::as_tibble(summary(sa_fit_EXclass_RE, pars = "class_mean"))
  cm_params <- cm_summary$parameter

  # Expected class means - shouldn't have sole-occupancy classes, except if these have shared SDs
  expected_cm <- setdiff(levels(sa_fit_EXclass_RE$network$classes),
                         c("Waitlist", "Pill placebo", "Psychological placebo"))

  expect_equal(cm_params, paste0("class_mean[", expected_cm, "]"))

  # Extract summaries for class_sd
  cs_summary <- tibble::as_tibble(summary(sa_fit_EXclass_RE, pars = "class_sd"))
  cs_params <- cs_summary$parameter

  # Expected class SDs - shared and named correctly, no sole-occupancy classes
  expected_cs <- forcats::fct_collapse(sa_fit_EXclass_RE$network$classes,
                                       `Exercise and SH no support` =
                                         c("Exercise promotion", "Self-help no support"),
                                       `SSRIs and NSSA` =
                                         c("SSRI/SNRI", "NSSA"),
                                       `Psychodynamic & Other psychological therapies` =
                                         c("Psychodynamic psychotherapy", "Other psychological therapies")) %>%
    levels() %>% setdiff(c("Waitlist", "Pill placebo", "Psychological placebo"))
  expect_equal(cs_params, paste0("class_sd[", expected_cs, "]"))
})
