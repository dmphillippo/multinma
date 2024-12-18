# Social Anxiety

library(multinma)
library(dplyr)
library(tidyr)
library(readr)

options(mc.cores = parallel::detectCores())

# Read in data
sa <- read_tsv("./data-raw/smoking/social_anxiety.txt")
sa_trt <- read_tsv("./data-raw/smoking/social_anxiety_treatments.txt")
sa_class <- read_tsv("./data-raw/smoking/social_anxiety_classes.txt")

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

social_anxiety <- as.data.frame(social_anxiety)

usethis::use_data(social_anxiety)

