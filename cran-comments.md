This release ensures compatibility with the forthcoming update to rstan and
StanHeaders, as well as providing extensive new features, improvements, and
fixes.

## Test environments
* local R installation (Windows 10), R 4.3.2
* win-builder (release, devel)
* Ubuntu 22.04.3 on GitHub Actions (release, devel, oldrel)
* Mac OS 12.7.2 on GitHub Actions (release)
* Windows Server 2022 10.0.20348 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 23 non-default packages.

* checking installed package size ... NOTE
  installed size is 18.2Mb
  sub-directories of 1Mb or more:
    R      1.6Mb
    doc    4.3Mb
    libs  11.6Mb

The docs directory contains 13 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

