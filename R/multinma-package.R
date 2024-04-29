#' multinma: A Package for Network Meta-Analysis of Individual and Aggregate
#' Data in Stan
#'
#' @description
#'   \if{html}{\figure{logo.svg}{options: width="120" style="float: right;" alt="multinma logo"}}
#'   An R package for performing network meta-analysis and network meta-regression
#'   with aggregate data, individual patient data, or mixtures of both.
#'
#' @details Network meta-analysis (NMA) combines (aggregate) data from multiple
#'   studies on multiple treatments in order to produce consistent estimates of
#'   relative treatment effects between each pair of treatments in the network
#'   \insertCite{TSD2}{multinma}.
#'
#'   Network meta-regression (NMR) extends NMA to include covariates, allowing
#'   adjustment for differences in effect-modifying variables between studies
#'   \insertCite{TSD3}{multinma}. NMR is typically performed using aggregate
#'   data (AgD), which lacks power and is prone to ecological bias. NMR with
#'   individual patient data (IPD) is the gold standard, if data are available.
#'
#'   Multilevel network meta-regression (ML-NMR) allows IPD and AgD to be
#'   incorporated together in a network meta-regression
#'   \insertCite{methods_paper,Phillippo_thesis}{multinma}. As in IPD NMR, an
#'   individual-level regression model is defined. AgD studies are then fitted
#'   by integrating the individual-level model over the respective covariate
#'   distributions. This correctly links the two levels of the model (instead of
#'   "plugging in" mean covariate values), avoiding aggregation bias.
#'   Population-adjusted treatment effects \insertCite{TSD18}{multinma} can be
#'   produced for any study population in the network, or for an external target
#'   population.
#'
#'   Models are estimated in a Bayesian framework using Stan
#'   \insertCite{Carpenter2017}{multinma}. Quasi-Monte Carlo numerical
#'   integration based on Sobol' sequences is used for the integration in ML-NMR
#'   models, with a Gaussian copula to account for correlations between
#'   covariates \insertCite{methods_paper,Phillippo_thesis}{multinma}.
#'
#' @section Getting Started:
#'   A good place to start is with the package vignettes which walk through
#'   example analyses, see `vignette("vignette_overview")` for an overview.
#'   The series of NICE Technical Support Documents on evidence synthesis gives
#'   a detailed introduction to network meta-analysis:
#'
#'   \insertRef{TSD_evsynth}{multinma}
#'
#'   Multilevel network meta-regression is set out in the following methods paper:
#'
#'   \insertRef{methods_paper}{multinma}
#'
#' @name multinma-package
#' @aliases multinma
#' @useDynLib multinma, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom dplyr %>%
#' @importFrom rlang abort warn inform enquo .data :=
#' @importFrom rstan sampling
#' @importFrom Rdpack reprompt
#' @importFrom graphics plot pairs
#' @importFrom grDevices nclass.Sturges
#' @importFrom stats complete.cases sd median quantile model.frame model.matrix
#'   model.offset terms optim pbinom dbinom qbinom as.formula update.formula
#'   weighted.mean runif dunif plogis pnorm qlogis qnorm uniroot update var
#'   setNames
#' @importFrom utils packageVersion head
#' @importFrom RcppParallel CxxFlags
#' @importFrom rstantools rstan_config
#'
#' @references
#'   \insertAllCited{}
#'
"_PACKAGE"

# Stop R CMD check thinking . used in pipes is an undeclared global variable
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

# Reexport survival::Surv()
#' @importFrom survival Surv
#' @export
survival::Surv
