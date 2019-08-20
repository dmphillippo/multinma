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
#' individual patient data from three studies (UNCOVER-1, -2, and
#' -3) \insertCite{Griffiths2015,Gordon2016}{multinma} and aggregate data from
#' one study (FIXTURE) \insertCite{Langley2014}{multinma}, respectively. Outcomes
#' are binary success/failure to achieve 75% reduction in symptoms on the
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
