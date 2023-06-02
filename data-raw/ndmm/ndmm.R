## Prepare `ndmm` datasets

# Simulated data to recreate Leahy and Walsh 2019

library(readr)
library(dplyr)
library(tidyr)
library(estmeansd)
library(simsurv)


# Set seed
set.seed(54321)

# Read in covariate summaries, survival distribution estimates
study_covs <- read_csv("./data-raw/ndmm/covariate_summaries.csv")  # Extracted from study reports
surv_pars <- read_csv("./data-raw/ndmm/genf_pars.csv")  # Estimated from digitized survival curves
surv_maxt_pcens <- read_csv("./data-raw/ndmm/maxt_pcens.csv")  # Calculated from digitized survival curves

# Reconstruct distribution for age
study_covs <- study_covs %>%
  rowwise() %>%
  mutate(age_fit = list(bc.mean.sd(age_min, age_iqr_l, age_median, age_iqr_u, age_max, sample_size)))

# Simulate covariates
ipd_all <-
  study_covs %>%
  group_by(study, treatment) %>%
  transmute(age = list(sample(age_fit[[1]]$bc.norm.rvs, size = sample_size, replace = FALSE)),
            iss_stage3 = list(rbinom(sample_size, 1, iss_stage3)),
            response_cr_vgpr = list(rbinom(sample_size, 1, response_cr_vgpr)),
            male = list(rbinom(sample_size, 1, male))) %>%
  unnest(cols = c(age, iss_stage3, response_cr_vgpr,
                  male)) %>%
  # Add treatment and study factor variables
  mutate(studyf = factor(study, levels = c("Attal2012", "McCarthy2012", "Palumbo2014", "Jackson2019", "Morgan2012")),
         trtf = factor(treatment, levels = c("Pbo", "Len", "Thal")),
         .after = treatment)

# Simulate survival outcomes using a generalised F model
Hgenf_reg <- function(t, x, betas, intercept, trteff, sigma, Q, P) {
  xtemp <- as.matrix(x[, names(betas)])
  flexsurv::Hgenf(t, mu = intercept + trteff + xtemp %*% betas, sigma = sigma, Q = Q, P = P)
}

# Prognostic and EM coefficients
betas <- -c(
  age = 0.065,
  iss_stage3 = 0.27,
  response_cr_vgpr = -0.2,
  male = 0.05,
  `age:trt` = -0.024,
  `iss_stage3:trt` = 0.26,
  `response_cr_vgpr:trt` = 0.31,
  `male:trt` = 0
)

surv_all <- ipd_all %>%
  group_by(study, treatment) %>%
  mutate(across(age:male, ~.x - mean(.x))) %>%
  mutate(`age:trt` = if_else(treatment != "Pbo", age, 0),
         `iss_stage3:trt` = if_else(treatment != "Pbo", iss_stage3, 0),
         `response_cr_vgpr:trt` = if_else(treatment != "Pbo", response_cr_vgpr, 0),
         `male:trt` = if_else(treatment != "Pbo", male, 0)) %>%
  nest(x = age:`male:trt`) %>%
  left_join(surv_pars) %>%
  left_join(surv_maxt_pcens) %>%
  mutate(trteff = recode(treatment,
                         Pbo = 0,
                         Len = treatmentLen,
                         Thal = treatmentThal),
         surv = list(simsurv(cumhazard = Hgenf_reg,
                             x = x[[1]],
                             betas = betas,
                             intercept = mu,
                             sigma = sigma,
                             P = P,
                             Q = Q,
                             trteff = trteff,
                             interval = c(1e-08, 10000),
                             maxt = max_time,
         ))
  ) %>%
  unnest(cols = c(x, surv)) %>%
  mutate(# Apply uniform censoring
    status = if_else(runif(n()) <= 0.2, 0L, status)
  )

# Create aggregate data for Morgan2012 and Jackson2016
agd_covs <- ipd_all %>% filter(studyf %in% c("Morgan2012", "Jackson2019")) %>%
  group_by(study, studyf, treatment, trtf) %>%
  summarise(sample_size = n(),
            as_tibble(as.list(setNames(fivenum(age), c("age_min", "age_iqr_l", "age_median", "age_iqr_h", "age_max")))),
            # Use estmeansd to calculate mean and sd
            setNames(as_tibble(unclass(estmeansd::bc.mean.sd(age_min, age_iqr_l, age_median, age_iqr_h, age_max, n = sample_size))[c("est.mean", "est.sd")]),
                     c("age_mean", "age_sd")),
            iss_stage3 = mean(iss_stage3),
            response_cr_vgpr = mean(response_cr_vgpr),
            male = mean(male)) %>%
  ungroup() %>%
  rename(trt = treatment)

agd_surv <- surv_all %>% filter(studyf %in% c("Morgan2012", "Jackson2019")) %>%
  ungroup() %>%
  transmute(study, studyf, trt = treatment, trtf, eventtime, status)


# Data for for IPD studies
ipd_surv <- ipd_all %>% filter(!studyf %in% c("Morgan2012", "Jackson2019")) %>%
  bind_cols(
    surv_all %>%
      ungroup() %>%
      filter(!studyf %in% c("Morgan2012", "Jackson2019")) %>%
      select(eventtime, status)
  ) %>%
  rename(trt = treatment) %>%
  ungroup()

# Name final datasets and store as standard data.frame
ndmm_ipd <- data.frame(ipd_surv)
ndmm_agd <- data.frame(agd_surv)
ndmm_agd_covs <- data.frame(agd_covs)

# Write datasets
usethis::use_data(ndmm_ipd, overwrite = TRUE)
usethis::use_data(ndmm_agd, overwrite = TRUE)
usethis::use_data(ndmm_agd_covs, overwrite = TRUE)
