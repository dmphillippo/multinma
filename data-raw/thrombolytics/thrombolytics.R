## Prepare `thrombolytics` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
thrombo_wide <- read_delim("./data-raw/thrombolytics/thrombolytics.txt", delim = " ")

# Make into long (tidy) format
thrombo <- thrombo_wide %>%
  group_by(study) %>%
  gather(key, var, r1:t3) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn = study,
            trtn = t,
            trtc = recode(trtn,
                          "SK",
                          "t-PA",
                          "Acc t-PA",
                          "SK + t-PA",
                          "r-PA",
                          "TNK",
                          "PTCA",
                          "UK",
                          "ASPAC"),
            r, n) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
thrombolytics <- as.data.frame(thrombo)

# Output to data dir ------------------------------------------------------
usethis::use_data(thrombolytics)
