Require StanHeaders version 2.32.9 or later in LinkingTo, as requested, which 
fixes the UBSAN warnings.

Confirmed that this fixes the sanitizer warnings with rocker/r-devel-san and 
rhub clang-asan.

## Test environments
* rocker/r-devel-san
* rhub clang-asan
* local R installation (Windows 10), R 4.4.0
* win-builder (release, devel)
* Ubuntu 22.04.4 on GitHub Actions (release, devel, oldrel)
* Mac OS 14.4.1 on GitHub Actions (release)
* Windows Server 2022 10.0.20348 on GitHub Actions (release)

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking package dependencies ... NOTE
  Imports includes 23 non-default packages.

* checking installed package size ... NOTE
  installed size is 18.0Mb
  sub-directories of 1Mb or more:
    R      1.7Mb
    doc    4.2Mb
    libs  11.5Mb
    
The docs directory contains 13 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.

## Downstream dependencies
There are currently no downstream dependencies for this package.

