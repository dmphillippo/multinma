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
#' Two data frames, `plaque_psoriasis_ipd` and `plaque_psoriasis_agd`,
#' containing (simulated) individual patient data from four studies and
#' aggregate data from five studies \insertCite{Phillippo_thesis}{multinma}.
#' Outcomes are binary success/failure to achieve 75%, 90%, or 100% reduction in
#' symptoms on the Psoriasis Area and Severity Index (PASI) scale.
#'
#' @format The individual patient data are contained in a data frame
#'   `plaque_psoriasis_ipd` with 4118 rows, one per individual, and 16 variables:
#'   \describe{
#'     \item{studyc}{study name}
#'     \item{trtc_long}{treatment name (long format)}
#'     \item{trtc}{treatment name}
#'     \item{trtn}{numeric treatment code}
#'     \item{pasi75}{binary PASI 75 outcome}
#'     \item{pasi90}{binary PASI 90 outcome}
#'     \item{pasi100}{binary PASI 100 outcome}
#'     \item{age}{age (years)}
#'     \item{bmi}{body mass index (BMI)}
#'     \item{pasi_w0}{PASI score at week 0}
#'     \item{male}{male sex (TRUE or FALSE)}
#'     \item{bsa}{body surface area (percent)}
#'     \item{weight}{weight (kilograms)}
#'     \item{durnpso}{duration of psoriasis (years)}
#'     \item{prevsys}{previous systemic treatment (TRUE or FALSE)}
#'     \item{psa}{psoriatic arthritis (TRUE or FALSE)}
#'   }
#'
#'  The aggregate data are contained in a data frame `plaque_psoriasis_agd` with 15
#'  rows, one per study arm, and 26 variables:
#'   \describe{
#'     \item{studyc}{study name}
#'     \item{trtc_long}{treatment name (long format)}
#'     \item{trtc}{treatment name}
#'     \item{trtn}{numeric treatment code}
#'     \item{pasi75_r, pasi75_n}{PASI 75 outcome count and denominator}
#'     \item{pasi90_r, pasi90_n}{PASI 75 outcome count and denominator}
#'     \item{pasi100_r, pasi100_n}{PASI 75 outcome count and denominator}
#'     \item{sample_size_w0}{sample size at week zero}
#'     \item{age_mean, age_sd}{mean and standard deviation of age (years)}
#'     \item{bmi_mean, bmi_sd}{mean and standard deviation of BMI}
#'     \item{pasi_w0_mean, pasi_w0_sd}{mean and standard deviation of PASI score at week 0}
#'     \item{male}{percentage of males}
#'     \item{bsa_mean, bsa_sd}{mean and standard deviation of body surface area (percent)}
#'     \item{weight_mean, weight_sd}{mean and standard deviation of weight (kilograms)}
#'     \item{durnpso_mean, durnpso_sd}{mean and standard deviation of duration of psoriasis (years)}
#'     \item{prevsys}{percentage of individuals with previous systemic treatment}
#'     \item{psa}{percentage of individuals with psoriatic arthritis}
#'   }
#'
#' @references \insertAllCited{}
#'
#' @rdname plaque_psoriasis
#' @aliases plaque_psoriasis
#'

"plaque_psoriasis_ipd"

#' @rdname plaque_psoriasis
"plaque_psoriasis_agd"

#' Thrombolytic treatments data
#'
#' Data frame containing the results of 50 trials of 8 thrombolytic drugs
#' (streptokinase, SK; alteplase, t-PA; accelerated alteplase, Acc t-PA;
#' streptokinase plus alteplase, SK+tPA; reteplase, r-PA; tenocteplase, TNK;
#' urokinase, UK; anistreptilase, ASPAC) plus per-cutaneous transluminal
#' coronary angioplasty (PTCA)
#' \insertCite{Boland2003,Lu2006,TSD4}{multinma}. The number of
#' deaths in 30 or 35 days following acute myocardial infarction are recorded.
#'
#' @format A data frame with 102 rows and 5 variables:
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

#' Beta blockers to prevent mortality after MI
#'
#' Data frame containing the number of deaths in 22 trials comparing beta
#' blockers vs. control for preventing mortality after myocardial infarction
#' \insertCite{Carlin1992,TSD2}{multinma}.
#'
#' @format A data frame with 44 rows and 5 variables:
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

"blocker"

#' Reduced dietary fat to prevent mortality
#'
#' Data frame containing the number of deaths and person-years at risk in 10
#' trials comparing reduced fat diets vs. control (non-reduced fat diet) for
#' preventing mortality \insertCite{Hooper2000,TSD2}{multinma}.
#'
#' @format A data frame with 21 rows and 7 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{studyc}{study name}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{r}{number of events}
#'   \item{n}{number randomised}
#'   \item{E}{person-years at risk}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"dietary_fat"

#' Mean off-time reduction in Parkison's disease
#'
#' Data frame containing the mean off-time reduction in patients given dopamine
#' agonists as adjunct therapy in Parkinson's disease, from 7 trials comparing
#' four active drugs and placebo \insertCite{TSD2}{multinma}.
#'
#' @format A data frame with 15 rows and 7 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{trtn}{numeric treatment code (placebo = 1)}
#'   \item{y}{mean off-time reduction}
#'   \item{se}{standard error}
#'   \item{n}{sample size}
#'   \item{diff}{mean difference vs. treatment in reference arm}
#'   \item{se_diff}{standard error of mean difference, see details}
#' }
#'
#' @details This dataset may be analysed using either an arm-based likelihood
#'   using `y` and `se`, or a contrast-based likelihood using `diff` and
#'   `se_diff` (or a combination of the two across different studies).
#'
#'   The contrast-based data is formatted as described in [set_agd_contrast()].
#'   That is, for the chosen reference arm in each study, the mean difference
#'   `diff` is set to `NA`, and `se_diff` is set to the standard error `se` of
#'   the outcome on the reference arm.
#'
#' @references
#'   \insertAllCited{}
#'

"parkinsons"

#' Stroke prevention in atrial fibrillation patients
#'
#' Data frame containing the results of 26 trials comparing 17 treatments in 4
#' classes for the prevention of stroke in patients with atrial fibrillation
#' \insertCite{Cooper2009}{multinma}. The data are the corrected versions
#' given by \insertCite{gemtc;textual}{multinma}.
#'
#' @format A data frame with 63 rows and 11 variables:
#' \describe{
#'   \item{studyc}{study name}
#'   \item{studyn}{numeric study ID}
#'   \item{trtc}{treatment name}
#'   \item{trtn}{numeric treatment code}
#'   \item{trt_class}{treatment class}
#'   \item{r}{number of events}
#'   \item{n}{sample size}
#'   \item{E}{person-years at risk}
#'   \item{stroke}{proportion of individuals with prior stroke}
#'   \item{year}{year of study publication}
#'   \item{followup}{mean length of follow-up (years)}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"atrial_fibrillation"

#' Statins for cholesterol lowering
#'
#' Data frame containing the results of 19 trials comparing statins to placebo
#' or usual care \insertCite{TSD3}{multinma}. The number of deaths (all-cause
#' mortality) are recorded. In some studies the aim was primary prevention
#' (patients had no previous heart disease), and in others the aim was secondary
#' prevention (patients had previous heart disease).
#'
#' @format A data frame with 38 rows and 7 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{studyc}{study name}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{prevention}{primary or secondary prevention study}
#'   \item{r}{number of deaths}
#'   \item{n}{sample size}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"statins"

#' BCG vaccination
#'
#' Data frame containing the results of 13 trials comparing BCG vaccination to
#' no vaccination for preventing tuberculosis (TB)
#' \insertCite{TSD3,Berkey1995}{multinma}. The numbers of individuals diagnosed
#' with TB in each arm during the study follow-up period are recorded. The
#' absolute degrees latitude at which the study was conducted are also recorded.
#'
#' @format A data frame with 26 rows and 6 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{latitude}{absolute degrees latitude}
#'   \item{r}{number diagnosed with TB}
#'   \item{n}{sample size}
#' }
#'
#' @references
#'   \insertAllCited{}
#'

"bcg_vaccine"

#' HTA Plaque Psoriasis
#'
#' Data frame containing the results of 16 trials comparing 8 treatments for
#' moderate-to-severe plaque psoriasis from an HTA report
#' \insertCite{Woolacott2006}{multinma}, analysed in TSD2
#' \insertCite{TSD2}{multinma}. Outcomes are success/failure to achieve 50%,
#' 75%, or 90% reduction in symptoms on the Psoriasis Area and Severity Index
#' (PASI) scale. Some studies report all three ordered outcomes, others only one
#' or two. The latter are coded as missing values (see details).
#'
#' @format A data frame with 36 rows and 9 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{studyc}{study name}
#'   \item{year}{year of publication}
#'   \item{trtn}{numeric treatment code}
#'   \item{trtc}{treatment name}
#'   \item{sample_size}{sample size in each arm}
#'   \item{PASI50, PASI75, PASI90}{ordered multinomial outcome counts (exclusive, see details)}
#' }
#'
#' @details Outcome counts are "exclusive"; that is, for a study reporting all
#'   outcomes, the counts represent the categories 50 < PASI < 75, 75 < PASI <
#'   90, and 90 < PASI < 100, and are named by the lower end of the interval.
#'   (As opposed to "inclusive" counts, which would represent the overlapping
#'   categories PASI > 50, PASI > 70, and PASI > 90.) The count for the fourth
#'   category (the lowest), 0 < PASI < 50, is equal to `sample_size - PASI50 -
#'   PASI75 - PASI90`.
#'
#'   Missing values are used where studies only report a subset of the outcomes.
#'   For a study reporting only two outcomes, say 50 and 75, the counts
#'   represent 50 < PASI < 75 and 75 < PASI < 100. For a study reporting only
#'   one outcome, say PASI 75, the count represents 75 < PASI < 100.
#'
#' @references
#'   \insertAllCited{}
#'

"hta_psoriasis"

#' Newly diagnosed multiple myeloma
#'
#' Three data frames, `ndmm_ipd`, `ndmm_agd`, and `ndmm_agd_covs` containing
#' (simulated) individual patient data (IPD) from three studies and aggregate
#' data (AgD) from two studies on newly diagnosed multiple myeloma. The outcome
#' of interest is progression-free survival after autologous stem cell
#' transplant. The IPD studies in `ndmm_ipd` provide event/censoring times and
#' covariate values for each individual. The AgD studies provide reconstructed
#' event/censoring times from digitized Kaplan-Meier curves in `ndmm_agd` and
#' covariate summaries in `ndmm_agd_covs`, obtained from published trial
#' reports. The data are constructed to resemble those used by
#' \insertCite{Leahy2019;textual}{multinma}.
#'
#' @format The individual patient data are contained in a data frame `ndmm_ipd`
#'   with 1325 rows, one per individual, and 10 variables:
#'   \describe{
#'     \item{study, studyf}{study name}
#'     \item{trt, trtf}{treatment name}
#'     \item{eventtime}{event/censoring time}
#'     \item{status}{censoring indicator (0 = censored, 1 = event)}
#'     \item{age}{age (years)}
#'     \item{iss_stage3}{ISS stage 3 (0 = no, 1 = yes)}
#'     \item{response_cr_vgpr}{complete or very good partial response (0 = no, 1 = yes)}
#'     \item{male}{male sex (0 = no, 1 = yes)}
#'   }
#'
#' @references \insertAllCited{}
#'
#' @rdname ndmm
"ndmm_ipd"

#' @rdname ndmm
#' @format The reconstructed Kaplan-Meier data for the aggregate studies are
#'   contained in a data frame `ndmm_agd` with 2819 rows and 6 variables:
#'  \describe{
#'    \item{study, studyf}{study name}
#'    \item{trt, trtf}{treatment name}
#'    \item{eventtime}{event/censoring time}
#'    \item{status}{censoring indicator (0 = censored, 1 = event)}
#'  }
"ndmm_agd"

#' @rdname ndmm
#' @format The covariate summaries extracted from published reportes for the
#'   aggregate studies are contained in a data frame `ndmm_agd_covs` with 4
#'   rows, one per study arm, and 15 columns:
#'  \describe{
#'    \item{study, studyf}{study name}
#'    \item{trt, trtf}{treatment name}
#'    \item{sample_size}{sample size in each arm}
#'    \item{age_min, age_iqr_l, age_median, age_iqr_h, age_max, age_mean, age_sd}{summary statistics for age (years)}
#'    \item{iss_stage3}{proportion of participants with ISS stage 3}
#'    \item{response_cr_vgpr}{proportion of participants with complete or very good partial response}
#'    \item{male}{proportion of male participants}
#'  }
"ndmm_agd_covs"

#' Social Anxiety
#'
#' Data frame containing the results of 101 clinical trials comparing 41 first-line treatments for social anxiety disorder
#' in adults . The 41 treatments are further categorised into 17 distinct classes. The interventions of interest include oral drugs,
#' psychological or behavioural therapies, and combinations of pharmacological and psychological therapies.
#' Following Mayo-Wilson et al. \insertCite{mayo2014psychological}{multinma} the data are given as standardised mean differences (SMD) to allow comparison across the different scales used in each study.
#'
#' @format A data frame with 248 rows and 8 variables:
#' \describe{
#'   \item{studyn}{numeric study ID}
#'   \item{studyc}{study name}
#'   \item{trtn}{treatment ID}
#'   \item{y}{standardised mean difference of the outcome measure}
#'   \item{se}{standard error of the outcome measure}
#'   \item{trtc}{treatment name}
#'   \item{classn}{class ID}
#'   \item{classc}{class name}
#' }
#'
#' @references
#'   \insertAllCited{}
#'
"social_anxiety"
