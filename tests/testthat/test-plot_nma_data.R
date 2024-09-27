
skip_on_cran()

library(vdiffr)
library(ggplot2)
library(dplyr)

test_that("Atrial fibrillation", {
  af_net <- set_agd_arm(atrial_fibrillation[atrial_fibrillation$studyc != "WASPO", ],
                        study = studyc,
                        trt = trtc,
                        r = r,
                        n = n,
                        trt_class = trt_class)

  expect_doppelganger("Atrial fibrillation network",
    plot(af_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE) +
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  )

  expect_doppelganger("Atrial fibrillation network star",
                      plot(af_net, layout = "star")
  )
})

test_that("BCG vaccine", {
  bcg_net <- set_agd_arm(bcg_vaccine,
                         study = studyn,
                         trt = trtc,
                         r = r,
                         n = n,
                         trt_ref = "Unvaccinated")

  expect_doppelganger("BCG vaccine network", plot(bcg_net))
})

test_that("Blocker", {
  blocker_net <- set_agd_arm(blocker,
                             study = studyn,
                             trt = trtc,
                             r = r,
                             n = n,
                             trt_ref = "Control")

  expect_doppelganger("Blocker network", plot(blocker_net))
})

test_that("Diabetes", {
  db_net <- set_agd_arm(diabetes,
                        study = studyc,
                        trt = trtc,
                        r = r,
                        n = n)
  expect_doppelganger("Diabetes network",
    plot(db_net, weight_edges = TRUE, weight_nodes = TRUE) + ggplot2::theme(legend.box.margin = ggplot2::unit(c(0, 0, 0, 4), "lines"))
  )
})

test_that("Dietary fat", {
  diet_net <- set_agd_arm(dietary_fat,
                          study = studyc,
                          trt = trtc,
                          r = r,
                          E = E,
                          trt_ref = "Control",
                          sample_size = n)

  expect_doppelganger("Dietary fat network", plot(diet_net))
})

test_that("HTA psoriasis", {
  pso_net <- set_agd_arm(hta_psoriasis,
                         study = paste(studyc, year),
                         trt = trtc,
                         r = multi(r0 = sample_size - rowSums(cbind(PASI50, PASI75, PASI90), na.rm = TRUE),
                                   PASI50, PASI75, PASI90,
                                   inclusive = FALSE,
                                   type = "ordered"))

  expect_doppelganger("HTA psoriasis network",
    plot(pso_net, weight_edges = TRUE, weight_nodes = TRUE) +
      # Nudge the legend over
      ggplot2::theme(legend.box.spacing = ggplot2::unit(0.75, "in"),
                     plot.margin = ggplot2::margin(0.1, 0, 0.1, 0.75, "in"))
  )
})

test_that("Parkinsons", {

  arm_net <- set_agd_arm(parkinsons,
                         study = studyn,
                         trt = trtn,
                         y = y,
                         se = se,
                         sample_size = n)

  expect_doppelganger("Parkinsons arm network",
    plot(arm_net, weight_edges = TRUE, weight_nodes = TRUE)
  )

  contr_net <- set_agd_contrast(parkinsons,
                                study = studyn,
                                trt = trtn,
                                y = diff,
                                se = se_diff,
                                sample_size = n)

  expect_doppelganger("Parkinsons contrast network",
    plot(contr_net, weight_edges = TRUE, weight_nodes = TRUE)
  )

  studies <- parkinsons$studyn
  (parkinsons_arm <- parkinsons[studies %in% 1:3, ])
  (parkinsons_contr <- parkinsons[studies %in% 4:7, ])

  mix_arm_net <- set_agd_arm(parkinsons_arm,
                             study = studyn,
                             trt = trtn,
                             y = y,
                             se = se,
                             sample_size = n)

  mix_contr_net <- set_agd_contrast(parkinsons_contr,
                                    study = studyn,
                                    trt = trtn,
                                    y = diff,
                                    se = se_diff,
                                    sample_size = n)

  mix_net <- combine_network(mix_arm_net, mix_contr_net)

  expect_doppelganger("Parkinsons mixed network",
    plot(mix_net, weight_edges = TRUE, weight_nodes = TRUE)
  )
})

test_that("Plaque psoriasis", {

  # IPD studies
  pso_ipd <- plaque_psoriasis_ipd %>%
    mutate(
      # Variable transformations
      bsa = bsa / 100,
      weight = weight / 10,
      durnpso = durnpso / 10,
      prevsys = as.numeric(prevsys),
      psa = as.numeric(psa),
      # Treatment classes
      trtclass = case_when(trtn == 1 ~ "Placebo",
                           trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                           trtn == 4 ~ "TNFa blocker",
                           trtn == 7 ~ "IL-12/23 blocker"),
      # Check complete cases for covariates of interest
      is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
    ) %>%
    arrange(studyc, trtn) %>%
    filter(is_complete)

  # AgD studies
  pso_agd <- plaque_psoriasis_agd %>%
    mutate(
      # Variable transformations
      bsa_mean = bsa_mean / 100,
      bsa_sd = bsa_sd / 100,
      weight_mean = weight_mean / 10,
      weight_sd = weight_sd / 10,
      durnpso_mean = durnpso_mean / 10,
      durnpso_sd = durnpso_sd / 10,
      prevsys = prevsys / 100,
      psa = psa / 100,
      # Treatment classes
      trtclass = case_when(trtn == 1 ~ "Placebo",
                           trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                           trtn == 4 ~ "TNFa blocker",
                           trtn == 7 ~ "IL-12/23 blocker")
    ) %>%
    arrange(studyc, trtn)

  pso_net <- combine_network(
    set_ipd(pso_ipd,
            study = studyc,
            trt = trtc,
            r = multi(r0 = 1,
                      PASI75 = pasi75,
                      PASI90 = pasi90,
                      PASI100 = pasi100,
                      type = "ordered", inclusive = TRUE),
            trt_class = trtclass),
    set_agd_arm(pso_agd,
                study = studyc,
                trt = trtc,
                r = multi(r0 = pasi75_n,
                          PASI75 = pasi75_r,
                          PASI90 = pasi90_r,
                          PASI100 = pasi100_r,
                          type = "ordered", inclusive = TRUE),
                trt_class = trtclass)
  )


  class_pal <- c("#D95F02", "#7570B3", "#E7298A", "#E6AB02")


  expect_doppelganger("Plaque psoriasis network",
    plot(pso_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE, nudge = 0.1) +
      ggraph::scale_edge_colour_manual("Data",
                                       values = c(AgD = "#113259", IPD = "#55A480")) +
      scale_fill_manual("Treatment class",
                        values = class_pal,
                        aesthetics = c("fill", "colour")) +
      guides(edge_colour = guide_legend(override.aes = list(edge_width = 2)),
             fill = guide_legend(override.aes = list(size = 2)))
  )
})

test_that("Smoking", {

  smknet <- set_agd_arm(smoking,
                        study = studyn,
                        trt = trtc,
                        r = r,
                        n = n,
                        trt_ref = "No intervention")

  expect_doppelganger("Smoking network",
    plot(smknet, weight_edges = TRUE, weight_nodes = TRUE) + ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 6), "lines"))
  )
})

test_that("Statins", {

  statin_net <- set_agd_arm(statins,
                            study = studyc,
                            trt = trtc,
                            r = r,
                            n = n,
                            trt_ref = "Placebo")

  expect_doppelganger("Statins network",
    plot(statin_net)
  )
})

test_that("Thrombolytics", {

  thrombo_net <- set_agd_arm(thrombolytics,
                             study = studyn,
                             trt = trtc,
                             r = r,
                             n = n)

  expect_doppelganger("Thrombolytics network",
    plot(thrombo_net, weight_edges = TRUE, weight_nodes = TRUE)
  )

  # Make trtf factor to order treatments in same way as Dias analysis
  trts <- dplyr::distinct(thrombolytics, trtn, trtc)
  trts <- dplyr::arrange(trts, trtn)
  thrombolytics$trtf <- factor(thrombolytics$trtn, levels = trts$trtn, labels = trts$trtc)
  thrombo_net2 <- set_agd_arm(thrombolytics,
                              study = studyn,
                              trt = trtf,
                              r = r,
                              n = n)

  expect_doppelganger("Thrombolytics Dias network",
                      plot(thrombo_net2, weight_edges = TRUE, weight_nodes = TRUE)
  )
})

test_that("Transfusion", {

  tr_net <- set_agd_arm(transfusion,
                        study = studyc,
                        trt = trtc,
                        r = r,
                        n = n,
                        trt_ref = "Control")

  expect_doppelganger("Transfusion network",
    plot(tr_net)
  )
})

test_that("NDMM", {

  ndmm_ipd$trtclass <- case_match(ndmm_ipd$trtf,
                                  "Pbo" ~ "Placebo",
                                  c("Len", "Thal") ~ "Active")

  ndmm_agd$trtclass <- case_match(ndmm_agd$trtf,
                                  "Pbo" ~ "Placebo",
                                  c("Len", "Thal") ~ "Active")

  ndmm_net <- combine_network(
    set_ipd(ndmm_ipd,
            study = studyf,
            trt = trtf,
            trt_class = trtclass,
            Surv = Surv(eventtime, status)),
    set_agd_surv(ndmm_agd,
                 study = studyf,
                 trt = trtf,
                 trt_class = trtclass,
                 Surv = Surv(eventtime, status),
                 covariates = ndmm_agd_covs)
  )

  ndmm_net <- add_integration(ndmm_net,
                              age = distr(qgamma, mean = age_mean, sd = age_sd),
                              iss_stage3 = distr(qbern, iss_stage3),
                              response_cr_vgpr = distr(qbern, response_cr_vgpr),
                              male = distr(qbern, male))

  expect_doppelganger("NDMM network",
    plot(ndmm_net,
         weight_nodes = TRUE,
         weight_edges = TRUE,
         # Nudge treatment labels away from nodes
         nudge = 0.1,
         # Manual layout
         layout = data.frame(x = c(0, -1, 1),
                             y = c(-0.5, 0, 0))) +
      guides(edge_colour = guide_legend(override.aes = list(edge_width = 2))) +
      theme(legend.position = "bottom", legend.direction = "vertical")
  )
})

test_that("Duplicated study arms plotted correctly", {

  dat <- data.frame(study = "a",
                    trt = c(1, 2, 2, 2),
                    r = 1, n = 1)

  net1 <- set_agd_arm(dat, study, trt, r = r, n = n)
  expect_doppelganger("Duplicated AgD arms network", plot(net1))


  net2 <- set_ipd(dat, study, trt, r = r)
  expect_doppelganger("Duplicated IPD arms network", plot(net2))

  dat_c <- data.frame(study = "a",
                      trt = c(1, 2, 2, 2),
                      y = c(NA, 1, 1, 1),
                      se = c(0.5, 1, 1, 1))
  net3 <- set_agd_contrast(dat_c, study, trt, y = y, se = se)
  expect_doppelganger("Duplicated AgD contrasts network", plot(net3))
})
