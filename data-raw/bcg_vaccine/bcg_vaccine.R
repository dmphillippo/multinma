## Prepare `bcg_vaccine` dataset

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
bcg_vaccine_wide <- read_tsv("./data-raw/bcg_vaccine/bcg_vaccine.txt") %>% mutate(studyn = 1:n())

# Make into long (tidy) format
bcg_vaccine <- bcg_vaccine_wide %>%
  group_by(studyn) %>%
  gather(key, var, matches("[rnt][0-9]+")) %>%
  mutate(arm = gsub("[rnt]([0-9]+)", "\\1", key),
         key = gsub("([rnt])[0-9]+", "\\1", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  transmute(studyn,
            trtn = t,
            trtc = recode(trtn,
                          "Unvaccinated",
                          "Vaccinated"),
            latitude,
            r, n) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
bcg_vaccine <- as.data.frame(bcg_vaccine)

# Add to package data directory
usethis::use_data(bcg_vaccine)
