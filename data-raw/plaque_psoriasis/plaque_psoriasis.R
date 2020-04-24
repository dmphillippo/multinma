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
         SEC_300 = 6L)
}

make_trtc_long <- function(trtc) {
  recode(trtc,
         PBO = "Placebo",
         IXE_Q2W = "Ixekizumab Q2W",
         IXE_Q4W = "Ixekizumab Q4W",
         ETN = "Etanercept",
         SEC_150 = "Secukinumab 150 mg",
         SEC_300 = "Secukinumab 300 mg")
}


# Individual patient data -------------------------------------------------
psoriasis_ipd <-
  read_csv("./data-raw/psoriasis/dummy_ipd.csv") %>%
  # Add in treatment codes
  mutate(trtn = make_trtn(trtc),
         trtc_long = make_trtc_long(trtc),
  # Add in numeric study IDs
         studyn = recode(study,
                         `UNCOVER-1` = 1L,
                         `UNCOVER-2` = 2L,
                         `UNCOVER-3` = 3L)) %>%
  rename(studyc = study,
         pasi75 = pasi75_w12_nri_01,
         prevsys = prevsys_01,
         psa = psa_01) %>%
  select(studyc, studyn, trtc_long, trtc, trtn, everything())

# Check treatment and study variables
distinct(psoriasis_ipd, trtc, trtn, trtc_long)
distinct(psoriasis_ipd, studyc, studyn)

# Store as standard data.frame
psoriasis_ipd <- as.data.frame(psoriasis_ipd)


# Aggregate data ----------------------------------------------------------
psoriasis_agd <- read_csv("./data-raw/psoriasis/fixture_agd.csv") %>%
  mutate(studyn = 4L,
         trtn = make_trtn(trtc),
         trtc_long = make_trtc_long(trtc)) %>%
  rename(studyc = study) %>%
  select(studyc, studyn, trtc_long, trtc, trtn, everything())

# Check treatment and study variables
distinct(psoriasis_agd, trtc, trtn, trtc_long)
distinct(psoriasis_agd, studyc, studyn)

# Store as standard data.frame
psoriasis_agd <- as.data.frame(psoriasis_agd)


# Output to data dir ------------------------------------------------------
usethis::use_data(psoriasis_ipd)
usethis::use_data(psoriasis_agd)
