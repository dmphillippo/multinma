This update fixes runtime errors in donttest examples (caused by latest
StanHeaders), as per email from CRAN.

## Test environments
* local R installation (Windows 10), R 4.3.0
* win-builder (release, devel)
* Ubuntu 22.04.2 on GitHub Actions (release, devel, oldrel)
* Mac OS 12.6.5 on GitHub Actions (release)
* Windows Server 2022 10.0.20348 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 22 non-default packages.

* checking installed package size ... 
  installed size is 12.3Mb
  sub-directories of 1Mb or more:
    doc    4.0Mb
    libs   6.9Mb

The docs directory contains 12 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

