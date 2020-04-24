## Prepare `psoriasis` datasets

library(dplyr)
library(readr)


# Treatment codings -------------------------------------------------------
make_trtn <- function(trtc) {
  recode(trtc,
         PBO = 1L,
         IXE_Q2W = 2L,
         IXE_Q4W = 3L,
         ETN = 4L,
         SEC_150 = 5L,
         SEC_300 = 6L,
         UST = 7L)
}

make_trtc_long <- function(trtc) {
  recode(trtc,
         PBO = "Placebo",
         IXE_Q2W = "Ixekizumab Q2W",
         IXE_Q4W = "Ixekizumab Q4W",
         ETN = "Etanercept",
         SEC_150 = "Secukinumab 150 mg",
         SEC_300 = "Secukinumab 300 mg",
         UST = "Ustekinumab")
}

YNtoTF <- function(x) {case_when(x %in% c("y", "Y") ~ TRUE, x %in% c("n", "N") ~ FALSE)}

# Individual patient data -------------------------------------------------
plaque_psoriasis_ipd <-
  read_csv("./data-raw/plaque_psoriasis/dummy_ipd.csv") %>%
  # Add in treatment codes
  mutate(trtn = make_trtn(trtc),
         trtc_long = make_trtc_long(trtc),
         sex = case_when(sex == "M" ~ TRUE, sex == "F" ~ FALSE)) %>%
  rename(studyc = study,
         pasi75 = pasi75_w12_nri,
         pasi90 = pasi90_w12_nri,
         pasi100 = pasi100_w12_nri,
         male = sex) %>%
  mutate_at(vars(pasi75, pasi90, pasi100, prevsys, psa), YNtoTF) %>%
  mutate_at(vars(pasi75, pasi90, pasi100), as.numeric) %>%
  select(studyc, trtc_long, trtc, trtn, pasi75, pasi90, pasi100, everything())

# Check treatment and variables
distinct(plaque_psoriasis_ipd, trtc, trtn, trtc_long)

# Store as standard data.frame
plaque_psoriasis_ipd <- as.data.frame(plaque_psoriasis_ipd)


# Aggregate data ----------------------------------------------------------
plaque_psoriasis_agd <- read_csv("./data-raw/plaque_psoriasis/all_agd.csv") %>%
  mutate(trtn = make_trtn(trtc),
         trtc_long = make_trtc_long(trtc)) %>%
  rename(studyc = study, male = sex) %>%
  select(studyc, trtc_long, trtc, trtn,
         pasi75_r, pasi75_n,
         pasi90_r, pasi90_n,
         pasi100_r, pasi100_n,
         everything())

# Check treatment and study variables
distinct(plaque_psoriasis_agd, trtc, trtn, trtc_long)

# Store as standard data.frame
plaque_psoriasis_agd <- as.data.frame(plaque_psoriasis_agd)


# Output to data dir ------------------------------------------------------
usethis::use_data(plaque_psoriasis_ipd)
usethis::use_data(plaque_psoriasis_agd)
