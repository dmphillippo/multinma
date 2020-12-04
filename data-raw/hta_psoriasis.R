## Prepare `hta_psoriasis` dataset from TSD2

library(dplyr)
library(tidyr)
library(readr)

# Read wide data
pso_wide <- read_tsv("./data-raw/psoriasis/psoriasis.txt") %>%
  mutate(studyn = 1:n())

# Make into long (tidy) format
pso <-
  pso_wide %>%
  group_by(studyn) %>%
  gather(key, var, matches("^[rnt]_")) %>%
  mutate(arm = gsub("[rnt]_([0-9]+)(_[0-9]+)?", "\\1", key),
         key = gsub("([rnt])_[0-9]+(_[0-9]+)?", "\\1\\2", key)) %>%
  spread(key, var) %>%
  filter(!is.na(t)) %>%
  ungroup() %>%
  mutate_at(vars(starts_with("C_")), ~if_else(is.na(.), 0, .)) %>%
  transmute(studyn,
            studyc = Trial,
            year = Year,
            trtn = t,
            trtc = recode(trtn,
                          "Supportive care",
                          "Etanercept 25 mg",
                          "Etanercept 50 mg",
                          "Efalizumab",
                          "Ciclosporin",
                          "Fumaderm",
                          "Infliximab",
                          "Methotrexate"),
            sample_size = n_1,
            # BEWARE: this code is not fully general, it is shortened to work with only the cases in this dataset
            PASI50 = case_when(C_2 == 2 ~ r_2, TRUE ~ NA_real_),
            PASI75 = case_when(C_2 == 2 & C_3 == 3 ~ r_3, C_2 == 3 ~ r_2, TRUE ~ NA_real_),
            PASI90 = case_when(C_3 == 3 & C_4 == 4 ~ sample_size - rowSums(cbind(r_1, r_2, r_3), na.rm = TRUE))) %>%
  mutate_at(vars(trtn, PASI50, PASI75, PASI90, sample_size), as.integer) %>%
  arrange(studyn, trtn)

# Store as standard data.frame
hta_psoriasis <- as.data.frame(pso)

usethis::use_data(hta_psoriasis, overwrite = TRUE)
