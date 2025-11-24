# multinma: A Package for Network Meta-Analysis of Individual and Aggregate Data in Stan

An R package for performing network meta-analysis and network
meta-regression with aggregate data, individual patient data, or
mixtures of both.

## Details

Network meta-analysis (NMA) combines (aggregate) data from multiple
studies on multiple treatments in order to produce consistent estimates
of relative treatment effects between each pair of treatments in the
network (Dias et al. 2011) .

Network meta-regression (NMR) extends NMA to include covariates,
allowing adjustment for differences in effect-modifying variables
between studies (Dias et al. 2011) . NMR is typically performed using
aggregate data (AgD), which lacks power and is prone to ecological bias.
NMR with individual patient data (IPD) is the gold standard, if data are
available.

Multilevel network meta-regression (ML-NMR) allows IPD and AgD to be
incorporated together in a network meta-regression (Phillippo et al.
2020; Phillippo 2019) . As in IPD NMR, an individual-level regression
model is defined. AgD studies are then fitted by integrating the
individual-level model over the respective covariate distributions. This
correctly links the two levels of the model (instead of "plugging in"
mean covariate values), avoiding aggregation bias. Population-adjusted
treatment effects (Phillippo et al. 2016) can be produced for any study
population in the network, or for an external target population.

Models are estimated in a Bayesian framework using Stan (Carpenter et
al. 2017) . Quasi-Monte Carlo numerical integration based on Sobol'
sequences is used for the integration in ML-NMR models, with a Gaussian
copula to account for correlations between covariates (Phillippo et al.
2020; Phillippo 2019) .

## Getting Started

A good place to start is with the package vignettes which walk through
example analyses, see
[`vignette("vignette_overview")`](https://dmphillippo.github.io/multinma/dev/articles/vignette_overview.md)
for an overview. The series of NICE Technical Support Documents on
evidence synthesis gives a detailed introduction to network
meta-analysis:

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Reken S, Ades AE
(2011). “NICE DSU Technical Support Documents 1-7: Evidence Synthesis
for Decision Making.” National Institute for Health and Care Excellence.
<https://www.sheffield.ac.uk/nice-dsu>.

Multilevel network meta-regression is set out in the following methods
paper:

Phillippo DM, Dias S, Ades AE, Belger M, Brnabic A, Schacht A, Saure D,
Kadziola Z, Welton NJ (2020). “Multilevel Network Meta-Regression for
population-adjusted treatment comparisons.” *Journal of the Royal
Statistical Society: Series A (Statistics in Society)*, **183**(3),
1189–1210. [doi:10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579)
.

## References

Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M,
Brubaker M, Guo J, Li P, Riddell A (2017). “Stan: A Probabilistic
Programming Language.” *Journal of Statistical Software*, **76**(1).
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01) .  
  
Dias S, Sutton AJ, Welton NJ, Ades AE (2011). “NICE DSU Technical
Support Document 3: Heterogeneity: subgroups, meta-regression, bias and
bias-adjustment.” National Institute for Health and Care Excellence.
<https://www.sheffield.ac.uk/nice-dsu>.  
  
Dias S, Welton NJ, Sutton AJ, Ades AE (2011). “NICE DSU Technical
Support Document 2: A generalised linear modelling framework for
pair-wise and network meta-analysis of randomised controlled trials.”
National Institute for Health and Care Excellence.
<https://www.sheffield.ac.uk/nice-dsu>.  
  
Phillippo DM (2019). *Calibration of Treatment Effects in Network
Meta-Analysis using Individual Patient Data*. Ph.D. thesis, University
of Bristol. Available from <https://research-information.bris.ac.uk/>.  
  
Phillippo DM, Ades AE, Dias S, Palmer S, Abrams KR, Welton NJ (2016).
“NICE DSU Technical Support Document 18: Methods for population-adjusted
indirect comparisons in submission to NICE.” National Institute for
Health and Care Excellence. <https://www.sheffield.ac.uk/nice-dsu>.  
  
Phillippo DM, Dias S, Ades AE, Belger M, Brnabic A, Schacht A, Saure D,
Kadziola Z, Welton NJ (2020). “Multilevel Network Meta-Regression for
population-adjusted treatment comparisons.” *Journal of the Royal
Statistical Society: Series A (Statistics in Society)*, **183**(3),
1189–1210. [doi:10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579)
.

## See also

Useful links:

- <https://dmphillippo.github.io/multinma/>

- <https://github.com/dmphillippo/multinma>

- Report bugs at <https://github.com/dmphillippo/multinma/issues>

## Author

**Maintainer**: David M. Phillippo <david.phillippo@bristol.ac.uk>
([ORCID](https://orcid.org/0000-0003-2672-7841))

Other contributors:

- Samuel J. Perren <samuel.perren@bristol.ac.uk>
  ([ORCID](https://orcid.org/0009-0005-1921-6957)) \[contributor\]
