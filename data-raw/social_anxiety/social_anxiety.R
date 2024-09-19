# Social Anxiety

library(multinma)
library(dplyr)
library(tidyr)
library(readr)

options(mc.cores = parallel::detectCores())

# Read in data
sa <- read_tsv("social_anxiety.txt")
sa_trt <- read_tsv("social_anxiety_treatments.txt")
sa_class <- read_tsv("social_anxiety_classes.txt")

# Tidy data
social_anxiety <- sa %>%
  mutate(Var1 = V) %>%
  pivot_longer(cols = matches("[0-9]+$"),
               names_to = c(".value", "arm"),
               names_pattern = "^([a-zA-Z]+)([0-9]+)$") %>%
  filter(!is.na(t)) %>%
  transmute(studyn, studyc, trtn = t, y = y, se = sqrt(Var)) %>%
  # Add in treatment names
  left_join(sa_trt, by = "trtn") %>%
  # Add in class details
  left_join(sa_class, by = "classn")

# Create network
sa_net <- set_agd_contrast(social_anxiety,
                           studyc, trtc,
                           y = y, se = se,
                           trt_class = classc)

sa_net

plot(sa_net, show_trt_class = TRUE)


# Fit RE NMA *without* any treatment classes
sa_fit_RE <- nma(sa_net,
                 trt_effects = "random",
                 prior_trt = normal(0, 100),
                 prior_het = half_normal(5))

sa_fit_RE

plot_prior_posterior(sa_fit_RE)

