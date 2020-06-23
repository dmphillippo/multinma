## Prepare `diabetes` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
diabetes_wide <- read_delim("./data-raw/diabetes/diabetes.txt", delim = " ")

# Make into long (tidy) format
diabetes <- diabetes_wide %>%
  mutate(studyn = 1:n()) %>%
  group_by(Study) %>%
  gather(key, var, t1:n3) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn,
            studyc = Study,
            trtn = t,
            trtc = recode(trtn,
                          "Diuretic",
                          "Placebo",
                          "Beta Blocker",
                          "CCB",
                          "ACE Inhibitor",
                          "ARB"),
            r, n, time) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
diabetes <- as.data.frame(diabetes)

# Output to data dir ------------------------------------------------------
usethis::use_data(diabetes)
