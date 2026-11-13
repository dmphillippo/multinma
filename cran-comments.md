Fixed UBSAN warnings identified in CRAN additional tests. Confirmed that this
fixes the warnings with rocker/r-devel-ubsan-clang.

## Test environments
* rocker/r-devel-ubsan-clang 
* local R installation (Windows 10), R 4.5.3
* Ubuntu 24.04.4 on GitHub Actions (release, devel, oldrel)
* Mac OS 15.7.4 on GitHub Actions (release)
* Windows Server 2025 10.0.26100 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 23 non-default packages.

* checking installed package size ... INFO
  installed size is 21.0Mb
  sub-directories of 1Mb or more:
    R      1.9Mb
    doc    4.8Mb
    libs  13.6Mb
    
The docs directory contains 15 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

