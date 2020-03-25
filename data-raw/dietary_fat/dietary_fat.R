## Prepare `dietary_fat` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
df_wide <- read_tsv("./data-raw/dietary_fat/dietary_fat.txt") %>% mutate(studyn = 1:n())

# Make into long (tidy) format
dietary_fat <- df_wide %>%
  group_by(studyn) %>%
  gather(key, var, t1:n3) %>%
  mutate(arm = gsub("[rntE]([0-9]+)", "\\1", key),
         key = gsub("([rntE])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn,
            studyc,
            trtn = t,
            trtc = recode(trtn,
                          "Control",
                          "Reduced Fat"),
            r, n, E) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
dietary_fat <- as.data.frame(dietary_fat)

# Add to package data directory
usethis::use_data(dietary_fat)
