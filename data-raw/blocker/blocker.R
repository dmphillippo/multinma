## Prepare `blocker` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
blocker_wide <- read_delim("./data-raw/blocker/blocker.txt", delim = " ")

# Make into long (tidy) format
blocker <- blocker_wide %>%
  mutate(studyn = 1:n()) %>%
  group_by(Study) %>%
  gather(key, var, r1:n2) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key),
         trtn = as.numeric(arm)) %>%
  spread(key, var) %>%
  filter(!is.na(trtn)) %>%
  ungroup() %>%
  transmute(studyn,
            trtn,
            trtc = recode(trtn,
                          "Control",
                          "Beta Blocker"),
            r, n) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
blocker <- as.data.frame(blocker)

# Output to data dir ------------------------------------------------------
usethis::use_data(blocker)
