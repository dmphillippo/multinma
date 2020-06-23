
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multinma: Network Meta-Analysis of individual and aggregate data in Stan

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/multinma)](https://CRAN.R-project.org/package=multinma)
[![Build
Status](https://travis-ci.org/dmphillippo/multinma.svg?branch=master)](https://travis-ci.org/dmphillippo/multinma)
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

## References

<div id="refs" class="references">

<div id="ref-Carpenter2017">

Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben
Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li,
and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.”
*Journal of Statistical Software* 76 (1).
<https://doi.org/10.18637/jss.v076.i01>.

</div>

<div id="ref-Phillippo_thesis">

Phillippo, David Mark. 2019. “Calibration of Treatment Effects in
Network Meta-Analysis Using Individual Patient Data.” PhD thesis,
University of Bristol.

</div>

<div id="ref-methods_paper">

Phillippo, David M., Sofia Dias, A. E. Ades, Mark Belger, Alan Brnabic,
Alexander Schacht, Daniel Saure, Zbigniew Kadziola, and Nicky J. Welton.
2020. “Multilevel Network Meta-Regression for Population-Adjusted
Treatment Comparisons.” *Journal of the Royal Statistical Society:
Series A (Statistics in Society)* 183 (3): 1189–1210.
<https://doi.org/10.1111/rssa.12579>.

</div>

</div>
