This is a minor bugfix release. Also addresses new CRAN check note "Rd \link{} 
targets missing package anchors".

## Test environments
* local R installation (Windows 10), R 4.4.1
* Ubuntu 22.04.4 on GitHub Actions (release, devel, oldrel)
* Mac OS 14.6.1 on GitHub Actions (release)
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

