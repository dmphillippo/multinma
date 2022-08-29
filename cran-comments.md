This update fixes NOTEs due to HTML5 checks on the manual, as per email from
CRAN.

## Test environments
* local R installation (Windows 10), R 4.2.1
* win-builder (release, devel)
* Ubuntu 20.04.4 on GitHub Actions (release, devel, oldrel)
* Mac OS 11.6.8 on GitHub Actions (release)
* Windows Server 2022 10.0.20348 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 4 notes

* checking CRAN incoming feasibility ... NOTE
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1002/jrsm.1167
    From: man/get_nodesplits.Rd
    Status: 503
    Message: Service Unavailable
    
  <... similar 503 errors for other DOIs ...>

Found the following (possibly) invalid DOIs:
  DOI: 10.1111/rssa.12579
    From: DESCRIPTION
          inst/CITATION
    Status: Service Unavailable
    Message: 503
    
The following DOIs are valid and load in a browser, but give false-positive 503
errors:
10.1177/0272989x9801800110
10.1111/rssa.12579
10.1002/sim.4780110202
10.1002/sim.3594
10.1002/sim.4780140406
10.1002/jrsm.1167
10.1002/sim.3767
10.1177/0272989X221117162

* checking package dependencies ... NOTE
  Imports includes 22 non-default packages.

* checking installed package size ... NOTE
  installed size is 12.3Mb
  sub-directories of 1Mb or more:
    R      1.1Mb
    doc    4.0Mb
    libs   6.7Mb

The docs directory contains 12 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

