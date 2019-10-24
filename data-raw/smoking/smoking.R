## Prepare `smoking` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
smk_wide <- read_tsv("./data-raw/smoking/smoking.txt") %>% mutate(studyn = 1:n())

# Make into long (tidy) format
smk <- smk_wide %>%
  group_by(studyn) %>%
  gather(key, var, r1:t3) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn,
            trtn = t,
            trtc = recode(trtn,
                          "No intervention",
                          "Self-help",
                          "Individual counselling",
                          "Group counselling"),
            r, n) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
smoking <- as.data.frame(smk)

# Add to package data directory
usethis::use_data(smoking)
