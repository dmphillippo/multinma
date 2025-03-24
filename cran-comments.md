This release adds new features and fixes a few bugs.

## Test environments
* local R installation (Windows 10), R 4.4.3
* Ubuntu 24.04.2 on GitHub Actions (release, devel, oldrel)
* Mac OS 14.7.4 on GitHub Actions (release)
* Windows Server 2022 10.0.20348 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 23 non-default packages.

* checking installed package size ... NOTE
  installed size is 22.2Mb
  sub-directories of 1Mb or more:
    R      2.1Mb
    doc    4.5Mb
    libs  14.8Mb
    
The docs directory contains 14 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

