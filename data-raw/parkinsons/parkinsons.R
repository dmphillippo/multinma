## Prepare `parkinsons` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
park_wide <- read_tsv("./data-raw/parkinsons/parkinsons.txt") %>% mutate(studyn = 1:n())

# Make into long (tidy) format
parkinsons <- park_wide %>%
  group_by(studyn) %>%
  gather(key, var, t1:n3) %>%
  mutate(arm = gsub("(t|y|n|se)([0-9]+)", "\\2", key),
         key = gsub("(t|y|n|se)([0-9]+)", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  # Compute diff, se_diff
  mutate(diff = if_else(arm == 1, NA_real_, y - first(y)),
         se_diff = if_else(arm == 1, se, sqrt(se^2 + first(se)^2))) %>%
  # Round to precision seen in TSD 2
  mutate(diff = round(diff, 2),
         se_diff = round(se_diff, 3)) %>%
  ungroup() %>%
  transmute(studyn,
            trtn = t,
            y, se, n,
            diff, se_diff) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
parkinsons <- as.data.frame(parkinsons)

usethis::use_data(parkinsons)
