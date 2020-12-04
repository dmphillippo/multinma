## Test environments
* local R installation (Windows 10), R 4.0.3
* win-builder (release, devel)
* Ubuntu 20.04 on GitHub Actions (release, devel)
* Mac OS 10.15.7 on GitHub Actions (release)
* Windows Server 2019 10.0.17763 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 22 non-default packages.

* checking installed package size ... NOTE
  installed size is 12.5Mb
  sub-directories of 1Mb or more:
    doc    5.7Mb
    libs   5.6Mb

The docs directory contains 12 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

