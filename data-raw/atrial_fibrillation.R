## Prepare "atrial_fibrillation" dataset

library(dplyr)

classes <- list(
  "Control" = c("01"),
  "Anti-coagulant" = c("02", "03", "04", "09"),
  "Anti-platelet" = c("05", "06", "07", "08", "10", "11", "12", "16", "17"),
  "Mixed" = c("13", "14", "15")
)

atrial_fibrillation <-
  left_join(gemtc::atrialFibrillation$data.ab,
            gemtc::atrialFibrillation$studies,
            by = "study") %>%
  left_join(gemtc::atrialFibrillation$treatments,
            by = c(treatment = "id")) %>%
  transmute(studyc = study,
            studyn = group_indices(., study),
            trtc = description,
            trtn = as.numeric(treatment),
            trt_class = forcats::fct_collapse(treatment, !!! classes),
            r = responders,
            n = sampleSize,
            E = exposure,
            stroke, year, followup)


atrial_fibrillation <- as.data.frame(atrial_fibrillation)

usethis::use_data(atrial_fibrillation)
