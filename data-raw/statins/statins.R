## Prepare `statins` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
statins_wide <- read_tsv("./data-raw/statins/statins.txt") %>% mutate(studyn = 1:n())

# Make into long (tidy) format
statins <- statins_wide %>%
  group_by(studyn) %>%
  gather(key, var, matches("[rnt][0-9]+")) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn,
            studyc,
            trtn = t,
            trtc = recode(trtn,
                          "Placebo",
                          "Statin"),
            prevention = recode(x,
                                "0" = "Primary",
                                "1" = "Secondary"),
            r, n) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
statins <- as.data.frame(statins)

# Add to package data directory
usethis::use_data(statins)
