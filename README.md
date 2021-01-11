
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multinma: Network Meta-Analysis of individual and aggregate data in Stan <img src='man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/multinma)](https://CRAN.R-project.org/package=multinma)
[![R build
status](https://github.com/dmphillippo/multinma/workflows/R-CMD-check/badge.svg)](https://github.com/dmphillippo/multinma/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3904454.svg)](https://doi.org/10.5281/zenodo.3904454)
<!-- badges: end -->

The `multinma` package implements network meta-analysis, network
meta-regression, and multilevel network meta-regression models which
combine evidence from a network of studies and treatments using either
aggregate data or individual patient data from each study (Phillippo et
al. 2020; Phillippo 2019). Models are estimated in a Bayesian framework
using Stan (Carpenter et al. 2017).

## Installation

You can install the released version of `multinma` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("multinma")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dmphillippo/multinma")
```

Installing from source (either from CRAN or GitHub) requires that the
`rstan` package is installed and configured. See the installation guide
[here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Getting started

A good place to start is with the package vignettes which walk through
example analyses, see `vignette("vignette_overview")` for an overview.
The series of NICE Technical Support Documents on evidence synthesis
gives a detailed introduction to network meta-analysis:

> Dias, S. et al. (2011). “NICE DSU Technical Support Documents 1-7:
> Evidence Synthesis for Decision Making.” *National Institute for
> Health and Care Excellence.* Available from <http://nicedsu.org.uk/>.

Multilevel network meta-regression is set out in the following methods
paper:

> Phillippo, D. M. et al. (2020). “Multilevel Network Meta-Regression
> for population-adjusted treatment comparisons.” *Journal of the Royal
> Statistical Society: Series A (Statistics in Society)*,
> 183(3):1189-1210. doi:
> [10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579).

## Citing multinma

The `multinma` package can be cited as follows:

> Phillippo, D. M. (2021). *multinma: Network Meta-Analysis of
> Individual and Aggregate Data in Stan*. R package version 0.2.1, doi:
> [10.5281/zenodo.3904454](https://doi.org/10.5281/zenodo.3904454).

When fitting ML-NMR models, please cite the methods paper:

> Phillippo, D. M. et al. (2020). “Multilevel Network Meta-Regression
> for population-adjusted treatment comparisons.” *Journal of the Royal
> Statistical Society: Series A (Statistics in Society)*,
> 183(3):1189-1210. doi:
> [10.1111/rssa.12579](https://doi.org/10.1111/rssa.12579).

## References

<div id="refs" class="references">

<div id="ref-Carpenter2017">

Carpenter, B., A. Gelman, M. D. Hoffman, D. Lee, B. Goodrich, M.
Betancourt, M. Brubaker, J. Guo, P. Li, and A. Riddell. 2017. “Stan: A
Probabilistic Programming Language.” *Journal of Statistical Software*
76 (1). <https://doi.org/10.18637/jss.v076.i01>.

</div>

<div id="ref-Phillippo_thesis">

Phillippo, D. M. 2019. “Calibration of Treatment Effects in Network
Meta-Analysis Using Individual Patient Data.” PhD thesis, University of
Bristol.

</div>

<div id="ref-methods_paper">

Phillippo, D. M., S. Dias, A. E. Ades, M. Belger, A. Brnabic, A.
Schacht, D. Saure, Z. Kadziola, and N. J. Welton. 2020. “Multilevel
Network Meta-Regression for Population-Adjusted Treatment Comparisons.”
*Journal of the Royal Statistical Society: Series A (Statistics in
Society)* 183 (3): 1189–1210. <https://doi.org/10.1111/rssa.12579>.

</div>

</div>
