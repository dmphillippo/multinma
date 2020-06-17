## Test environments
* local R installation (windows), R 4.0.1, 4.0.0
* win-builder (devel, release)
* rhub (ubuntu release, fedora devel)

## R CMD check results

0 errors | 0 warnings | 4 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'David M. Phillippo <david.phillippo@bristol.ac.uk>'

New submission

This is the first submission of this package.

* checking package dependencies ... NOTE
Imports includes 25 non-default packages.
Importing from so many packages makes the package vulnerable to any of
them becoming unavailable.  Move as many as possible to Suggests and
use conditionally.

* checking installed package size ... NOTE
  installed size is 20.8Mb
  sub-directories of 1Mb or more:
    doc   10.1Mb
    libs   9.7Mb

The docs directory contains 11 html vignettes demonstrating full analyses using
the package.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

Required to compile Stan models with rstan package.
