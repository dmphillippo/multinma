# Overview of Examples

This package contains a number of vignettes, each one walking through an
example analysis. The table below gives an overview.

Many of these examples recreate analyses from the series of Technical
Support Documents published by the NICE Decision Support Unit ([Dias et
al. 2011](#ref-TSD_evsynth)). The exceptions are atrial fibrillation
([Cooper et al. 2009](#ref-Cooper2009)), white blood cell transfusion
([Turner et al. 2012](#ref-Turner2012)), social anxiety ([Perren et al.
2025](#ref-Perren2025); [Mayo-Wilson et al.
2014](#ref-mayo2014psychological)), and plaque psoriasis multilevel
network meta-regression ([Phillippo et al. 2020](#ref-methods_paper),
[2022](#ref-Phillippo2022)).

| Title                                                                                                          | Outcome type                                           | Likelihood                                                                | Link function | Notable features                                                              |
|:---------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------|:--------------------------------------------------------------------------|:--------------|:------------------------------------------------------------------------------|
| [Blocker](https://dmphillippo.github.io/multinma/dev/articles/example_blocker.md)                              | Counts                                                 | Binomial                                                                  | logit         | Pairwise MA                                                                   |
| [Dietary fat](https://dmphillippo.github.io/multinma/dev/articles/example_dietary_fat.md)                      | Rates                                                  | Poisson                                                                   | log           | Analysis of log rate ratios from rate data                                    |
| [Diabetes](https://dmphillippo.github.io/multinma/dev/articles/example_diabetes.md)                            | Counts with time at risk                               | Binomial                                                                  | cloglog       | Analysis of log hazard ratios from count data with time at risk               |
| [Parkinson’s](https://dmphillippo.github.io/multinma/dev/articles/example_parkinsons.md)                       | Continuous                                             | Normal                                                                    | Identity      | Analysis of arm-based data, contrast-based data, and a mixture of both        |
| [HTA plaque psoriasis](https://dmphillippo.github.io/multinma/dev/articles/example_hta_psoriasis.md)           | Ordered                                                | Multinomial (ordered)                                                     | probit        | Analysis of ordered categorical outcomes                                      |
| [Statins](https://dmphillippo.github.io/multinma/dev/articles/example_statins.md)                              | Counts                                                 | Binomial                                                                  | logit         | Meta-regression with subgroups                                                |
| [BCG vaccine](https://dmphillippo.github.io/multinma/dev/articles/example_bcg_vaccine.md)                      | Counts                                                 | Binomial                                                                  | logit         | Meta-regression with a continuous covariate, predictive distributions         |
| [Smoking cessation](https://dmphillippo.github.io/multinma/dev/articles/example_smoking.md)                    | Counts                                                 | Binomial                                                                  | logit         | Assessing inconsistency with unrelated mean effects and node-splitting models |
| [Thrombolytics](https://dmphillippo.github.io/multinma/dev/articles/example_thrombolytics.md)                  | Counts                                                 | Binomial                                                                  | logit         | Assessing inconsistency with unrelated mean effects and node-splitting models |
| [Atrial fibrillation](https://dmphillippo.github.io/multinma/dev/articles/example_atrial_fibrillation.md)      | Counts                                                 | Binomial                                                                  | logit         | Meta-regression with shared class interactions                                |
| [WBC transfusion](https://dmphillippo.github.io/multinma/dev/articles/example_transfusion.md)                  | Counts                                                 | Binomial                                                                  | logit         | Informative log-Normal prior on \tau^2                                        |
| [ML-NMR plaque psoriasis](https://dmphillippo.github.io/multinma/dev/articles/example_plaque_psoriasis.md)     | Binary (IPD) and counts (AgD), and ordered categorical | Bernoulli (IPD) and two-parameter Binomial (AgD), and ordered multinomial | probit        | Multilevel network meta-regression combining IPD and AgD                      |
| [ML-NMR newly diagnosed multiple myeloma](https://dmphillippo.github.io/multinma/dev/articles/example_ndmm.md) | Time-to-event with censoring                           | M-spline baseline hazard                                                  | log           | Multilevel network meta-regression combining IPD and AgD                      |
| [Social anxiety](https://dmphillippo.github.io/multinma/dev/articles/example_social_anxiety.md)                | Continuous                                             | Normal                                                                    | Identity      | Model selection with class effects models                                     |
| [Certolizumab](https://dmphillippo.github.io/multinma/dev/articles/example_certolizumab.md)                    | Counts                                                 | Binomial                                                                  | logit         | Baseline risk meta-regression                                                 |

## References

Cooper, N. J., A. J. Sutton, D. Morris, A. E. Ades, and N. J. Welton.
2009. “Addressing Between-Study Heterogeneity and Inconsistency in Mixed
Treatment Comparisons: Application to Stroke Prevention Treatments in
Individuals with Non-Rheumatic Atrial Fibrillation.” *Statistics in
Medicine* 28 (14): 1861–81. <https://doi.org/10.1002/sim.3594>.

Dias, S., N. J. Welton, A. J. Sutton, D. M. Caldwell, G. Lu, S. Reken,
and A. E. Ades. 2011. “NICE DSU Technical Support Documents 1-7:
Evidence Synthesis for Decision Making.” National Institute for Health
and Care Excellence. <https://sheffield.ac.uk/nice-dsu>.

Mayo-Wilson, Evan, Sofia Dias, Ifigeneia Mavranezouli, Kayleigh Kew,
David M Clark, AE Ades, and Stephen Pilling. 2014. “Psychological and
Pharmacological Interventions for Social Anxiety Disorder in Adults: A
Systematic Review and Network Meta-Analysis.” *The Lancet Psychiatry* 1
(5): 368–76.

Perren, Samuel J., Hugo Pedder, Nicky J. Welton, and David M. Phillippo.
2025. “Network Meta-Analysis with Class Effects: A Practical Guide and
Model Selection Algorithm.” *Medical Decision Making* 46 (3): 275–95.
<https://doi.org/10.1177/0272989x251389887>.

Phillippo, D. M., S. Dias, A. E. Ades, M. Belger, A. Brnabic, D. Saure,
Y. Schymura, and N. J. Welton. 2022. “Validating the Assumptions of
Population Adjustment: Application of Multilevel Network Meta-Regression
to a Network of Treatments for Plaque Psoriasis.” *Medical Decision
Making*. <https://doi.org/10.1177/0272989X221117162>.

Phillippo, D. M., S. Dias, A. E. Ades, M. Belger, A. Brnabic, A.
Schacht, D. Saure, Z. Kadziola, and N. J. Welton. 2020. “Multilevel
Network Meta-Regression for Population-Adjusted Treatment Comparisons.”
*Journal of the Royal Statistical Society: Series A (Statistics in
Society)* 183 (3): 1189–1210. <https://doi.org/10.1111/rssa.12579>.

Turner, R. M., J. Davey, M. J. Clarke, S. G. Thompson, and J. P. T.
Higgins. 2012. “Predicting the Extent of Heterogeneity in Meta-Analysis,
Using Empirical Data from the Cochrane Database of Systematic Reviews.”
*International Journal of Epidemiology* 41 (3): 818–27.
<https://doi.org/10.1093/ije/dys041>.
