Incremented the version number and resubmitted.

Package previously archived due to UBSAN warnings. StanHeaders 2.32.7 and 2.32.8
omitted the necessary patch, which has been reinstated with 2.32.9.

Confirmed that this fixes the sanitizer warnings with rocker/r-devel-san. No 
changes are required to this package, this submission is identical to the 
previous release (v0.7.0).

## Test environments
* rocker/r-devel-san
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

