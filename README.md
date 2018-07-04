# polyIBD
[![Travis-CI Build Status](https://travis-ci.org/nickbrazeau/polyIBD.svg?branch=master)](https://travis-ci.org/nickbrazeau/polyIBD)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nickbrazeau/polyIBD?branch=master&svg=true)](https://ci.appveyor.com/project/nickbrazeau/polyIBD)
[![Coverage Status](https://img.shields.io/codecov/c/github/nickbrazeau/polyIBD/master.svg)](https://codecov.io/github/nickbrazeau/polyIBD?branch=master)

## Status
This is the branch that was presented at the Genomic Epidemiology of Malaria (GEM) Conference on June 13, 2018. We have now branched to `devo` work on the within relatedness (F<sub>ws</sub>) parameter.  
  
For the time being, the `master` branch will be stable but is considered alpha-version of the software. 

## Purpose 
  
R package for inferring blocks of identity by descent (IBD) from un-phased haplotypic data.
  

## Installation 
To install the current version of _polyIBD_, first ensure that you have the R-package _devtools_:
``` r
install.packages("devtools")
```

Then download and build _polyIBD_ from Github with: 
``` r
devtools::install_github("nickbrazeau/polyIBD")
```
