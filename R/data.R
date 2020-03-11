#' Smoking cessation data
#'
#' Data frame containing the results of 24 trials of 4 smoking cessation
#' treatments \insertCite{Hasselblad1998,TSD4}{multinma}.
#'
#' @format A data frame with 50 rows and 5 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{r}{total number of events}
#'   \item{n}{total number of individuals}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"smoking"

#' Plaque psoriasis data
#'
#' Two data frames, `psoriasis_ipd` and `psoriasis_agd`, containing (simulated)
#' individual patient data from three studies (UNCOVER-1, -2, and -3)
#' \insertCite{Griffiths2015,Gordon2016}{multinma} and aggregate data from one
#' study (FIXTURE) \insertCite{Langley2014}{multinma}, respectively. Outcomes
#' are binary success/failure to achieve 75\% reduction in symptoms on the
#' Psoriasis Area and Severity Index (PASI) scale.
#'
#' @format The individual patient data are contained in a data frame
#'   `psoriasis_ipd` with 3858 rows, one per individual, and 11 variables:
#'   \describe{
#'     \item{studyc}{study name}
#'     \item{studyn}{numeric study ID}
#'     \item{trtc_long}{treatment name (long format)}
#'     \item{trtc}{treatment name}
#'     \item{trtn}{numeric treatment code}
#'     \item{bsa}{body surface area (percent)}
#'     \item{weight}{weight (kilograms)}
#'     \item{durnpso}{duration of psoriasis (years)}
#'     \item{prevsys}{previous systemic treatmet (yes = 1, no = 0)}
#'     \item{psa}{psoriatic arthritis (yes = 1, no = 0)}
#'     \item{pasi75}{binary PASI 75 outcome (success = 1, failure = 0)}
#'   }
#'
#'  The aggregate data are contained in a data frame `psoriasis_agd` with 4
#'  rows, one per study arm, and 16 variables:
#'   \describe{
#'     \item{studyc}{study name}
#'     \item{studyn}{numeric study ID}
#'     \item{trtc_long}{treatment name (long format)}
#'     \item{trtc}{treatment name}
#'     \item{trtn}{numeric treatment code}
#'     \item{sample_size_w0}{sample size at week zero}
#'     \item{bsa_mean}{mean body surface area (percent)}
#'     \item{bsa_sd}{standard deviation of body surface area (percent)}
#'     \item{weight_mean}{mean weight (kilograms)}
#'     \item{weight_sd}{standard deviation of weight (kilograms)}
#'     \item{durnpso_mean}{mean duration of psoriasis (years)}
#'     \item{durnpso_sd}{standard deviation of duration of psoriasis (years)}
#'     \item{prevsys}{percentage of individuals with previous systemic treatment}
#'     \item{psa}{percentage of individuals with psoriatic arthritis}
#'     \item{pasi75_n}{PASI 75 outcome denominator}
#'     \item{pasi75_r}{PASI 75 outcome numerator}
#'   }
#'
#' @references \insertAllCited{}
#'
#' @rdname psoriasis
#' @aliases psoriasis psoriasis_agd
#'

"psoriasis_ipd"

#' Thrombolytic treatments data
#'
#' Data frame containing the results of 50 trials of 8 thrombolytic drugs
#' (streptokinase, SK; alteplase, t-PA; accelerated alteplase, Acc t-PA;
#' streptokinase plus alteplase, SK+tPA; reteplase, r-PA; tenocteplase, TNK;
#' urokinase, UK; anistreptilase, ASPAC) plus per-cutaneous transluminal
#' coronary angioplasty (PTCA)
#' \insertCite{Boland2003,Lu2006,TSD4}{multinma}. The number of
#' deaths in 30 or 35 days following accute myocardial infarction are recorded.
#'
#' @format A data frame with 50 rows and 5 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{r}{total number of events}
#'   \item{n}{total number of individuals}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"thrombolytics"

#' Incidence of diabetes in trials of antihypertensive drugs
#'
#' Data frame containing the number of new cases of diabetes in 22 trials of 6
#' antihypertensive drugs \insertCite{Elliott2007,TSD2}{multinma}. The trial
#' duration (in years) is also recorded.
#'
#' @format A data frame with 48 rows and 7 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{studyc}{study name}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{r}{total number of events}
#'   \item{n}{total number of individuals}
#'   \item{time}{trial follow-up (years)}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"diabetes"

#' Granulocyte transfusion in patients with neutropenia or neutrophil
#' dysfunction
#'
#' Data frame containing the number of deaths in 6 trials comparing transfusion
#' of granulocytes (white blood cells) to control
#' \insertCite{Stanworth2005}{multinma}. Previously used to demonstrate
#' informative prior distributions for the heterogeneity variance by
#' \insertCite{Turner2012;textual}{multinma}.
#'
#' @format A data frame with 12 rows and 4 variables:
#' \describe{
#'   \item{studyc}{study name}
#'   \item{trtc}{treatment name}
#'   \item{r}{total number of deaths}
#'   \item{n}{total number of individuals}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"transfusion"
