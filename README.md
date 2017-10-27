[![Build Status](https://travis-ci.org/MahShaaban/pcr.svg?branch=master)](https://travis-ci.org/MahShaaban/pcr)
[![Build status](https://ci.appveyor.com/api/projects/status/y9hfiwwc390cce28?svg=true)](https://ci.appveyor.com/project/MahShaaban/pcr)
[![Build status](https://ci.appveyor.com/api/projects/status/y9hfiwwc390cce28/branch/master?svg=true)](https://ci.appveyor.com/project/MahShaaban/pcr/branch/master)
[![Coverage Status](https://img.shields.io/codecov/c/github/MahShaaban/pcr/master.svg)](https://codecov.io/github/MahShaaban/pcr?branch=master)

# pcr  
## Overview  

## What is pcr used for?  

## Getting started  


```r
# install package from github (under development)
devtools::install_github('MahShaaban/pcr')
```
```r
# load required libraries
library(pcr)
library(ggplot2)
```

```r
# default mode delta_delta_ct
## locate and read raw ct data
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- readr::read_csv(fl)

## add grouping variable
group_var <- rep(c('brain', 'kidney'), each = 6)

# calculate all values and errors in one step
## mode == 'average_ct' default
res <- pcr_analyze(ct1,
  group_var = group_var,
  reference_gene = 'GAPDH',
  reference_group = 'brain')
res
```

```r
## locate and read data
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- readr::read_csv(fl)

## make a vector of RNA amounts
amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

## calculate amplification effeciency
res <- pcr_assess(ct3,
                  amount = amount,
                  reference_gene = 'GAPDH',
                  method = 'effeciency')
res
```
### Other analysis methods  

* Delta CT method  
* Relative standard curve method  

### Testing statistical significance  
* Two group testing **t-**test and **wilcoxon** test  
* Linear regression testing  

## Documnetation  

```r
browseVignettes("pcr")
```  

## Citation  

```r
citation("pct")
```
