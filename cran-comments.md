This update fixes an invalid img tag height attribute in the documentation, as
per email from CRAN.

## Test environments
* local R installation (Windows 10), R 4.1.2
* win-builder (release, devel)
* Ubuntu 20.04.3 on GitHub Actions (release, devel, oldrel)
* Mac OS 11.6.2 on GitHub Actions (release)
* Windows Server 2019 10.0.17763 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 22 non-default packages.

* checking installed package size ... NOTE
  installed size is 9.7Mb
  sub-directories of 1Mb or more:
    doc    3.7Mb
    libs   6.0Mb

The docs directory contains 12 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

