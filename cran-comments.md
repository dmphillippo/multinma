This release addresses an "additional checks" UBSAN failure, caused by a memory
allocation bug in Stan (https://github.com/stan-dev/rstan/issues/1111).
StanHeaders has now been patched which resolves the issue; no changes are
necessary in multinma. Verified passing UBSAN checks with rocker/r-devel-san.

In response to CRAN request, I have also:

* Added the few missing \value tags to documentation

Regarding the other issues raised:

* \dontrun{} is used in man/pairs.stan_nma.Rd, as the example is designed to
  raise warning messages. These are divergent transitions in Stan, which we then
  demonstrate for the user how to examine using the pairs() plot.
  
* The global environment is never modified. rm(list = ls()) in some large test
  files is never run on CRAN or by the user, these lines only run when testing 
  on GitHub Actions to avoid running out of memory. These tests are skipped on 
  CRAN to keep check time manageable.
  
* Similarly, .GlobalEnv is never modified. In the .Rd-files starting with 
  aa_example_, the assign() to .GlobalEnv only runs and is only needed when 
  building the package website using pkgdown, i.e. when pkgdown::in_pkgdown() is
  TRUE. This code is never run otherwise.
  

## Test environments
* local rocker/r-devel-san image with sanitizers
* local R installation (Windows 10), R 4.3.3
* win-builder (release, devel)
* Ubuntu 22.04.4 on GitHub Actions (release, devel, oldrel)
* Mac OS 12.7.3 on GitHub Actions (release)
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

