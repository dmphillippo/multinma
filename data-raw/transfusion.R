## Prepare "transfusion" dataset

library(dplyr)

transfusion <- tribble(
  ~studyc, ~trtc, ~r, ~n,
  "Bow 1984", "Transfusion", 5, 13,
  "Bow 1984", "Control", 4, 11,
  "Herzig 1977", "Transfusion", 1, 13,
  "Herzig 1977", "Control", 3, 14,
  "Higby 1975", "Transfusion", 2, 17,
  "Higby 1975", "Control", 14, 19,
  "Scali 1978", "Transfusion", 0, 13,
  "Scali 1978", "Control", 1, 12,
  "Vogler 1977", "Transfusion", 7, 17,
  "Vogler 1977", "Control", 9, 13,
  "Winston 1982a", "Transfusion", 18, 48,
  "Winston 1982a", "Control", 13, 47,
)

transfusion <- as.data.frame(transfusion)

usethis::use_data(transfusion)
