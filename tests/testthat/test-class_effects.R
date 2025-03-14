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
                                          prior_class_sd = normal(0.33,0.1),
                                          class_sd =
                                            list(`Exercise and SH no support` = c("Exercise promotion", "Self-help no support"),
                                                 `SSRIs and NSSA` = c("SSRI/SNRI", "NSSA"),
                                                 `Psychodynamic & Other psychological therapies` = c("Psychodynamic psychotherapy", "Other psychological therapies")),
                                          test_grad = TRUE))

test_that("class_mean and class_sd parameters are named correctly", {
  # Extract summaries for class_mean
  cm_summary <- tibble::as_tibble(summary(sa_fit_EXclass_RE, pars = "class_mean"))
  cm_params <- cm_summary$parameter

  expected_possible_cm <- paste0("class_mean[", unique(sa_fit_EXclass_RE$network$classes), "]")

  #Check that every actual parameter name is one of the expected possible names
  expect_true(all(cm_params %in% expected_possible_cm),
              info = "Each class_mean parameter should be constructible from the unique classes in the network.")

  # Extract summaries for class_sd
  cs_summary <- tibble::as_tibble(summary(sa_fit_EXclass_RE, pars = "class_sd"))
  cs_params <- cs_summary$parameter

  expected_possible_cs <- paste0("class_sd[", unique(sa_fit_EXclass_RE$network$class_sd), "]")

  #Check that every actual parameter name is one of the expected possible names
  expect_true(all(cs_params %in% expected_possible_cs),
              info = "Each class_sd parameter should be constructible from the class_sd list in the network.")
})
