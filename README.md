[![Build Status](https://travis-ci.org/MahShaaban/pcr.svg?branch=master)](https://travis-ci.org/MahShaaban/pcr)
[![Build status](https://ci.appveyor.com/api/projects/status/y9hfiwwc390cce28?svg=true)](https://ci.appveyor.com/project/MahShaaban/pcr)
[![Build status](https://ci.appveyor.com/api/projects/status/y9hfiwwc390cce28/branch/master?svg=true)](https://ci.appveyor.com/project/MahShaaban/pcr/branch/master)
[![Coverage Status](https://img.shields.io/codecov/c/github/MahShaaban/pcr/master.svg)](https://codecov.io/github/MahShaaban/pcr?branch=master)
[![CRAN version](https://img.shields.io/badge/CRAN-v1.1.2-blue.svg)](https://CRAN.R-project.org/package=pcr) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/pcr)  

# pcr  

# Overview  

Quantitative real-time PCR is an imprtant technique in medical and biomedical applicaitons. The `pcr` package provides a unified interface for quality assessing, analyzing and testing qPCR data for statistical significance. The aim of this document is to describe the different methods and modes used to relatively quantify gene expression of qPCR and their implemenation in the `pcr` package.  

# Getting started 

The `pcr` is available on CRAN. To install it, use:  

```r
# install package CRAN
install.packages('pcr')
```

The development version of the package can be obtained through:  

```r
# install package from github (under development)
devtools::install_github('MahShaaban/pcr')
```

```r
# load required libraries
library(pcr)
```

The following chunck of code locates a dataset of CT values of two genes from 12 different samples and performs a quick analysis to obtain the expression of a target gene **c-myc** normalized by a control GAPDH in the **Kidney** samples relative to the brain samples. `pcr_analyze` provides differnt methods, the default one that is used here is 'delta_delta_ct' applies the popular Double Delta CT method.  

```r
# default mode delta_delta_ct
## locate and read raw ct data
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- readr::read_csv(fl)

## add grouping variable
group_var <- rep(c('brain', 'kidney'), each = 6)

# calculate all values and errors in one step
## mode == 'separate_tube' default
res <- pcr_analyze(ct1,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain')
  
res
```

The output of `pcr_analyze` is explained in the documnetation of the function `?pcr_analyze` and the method it calls `?pcr_ddct`. Briefly, the input includes the CT value of c-myc `normalized` to the control GAPDH, The `calibrated` value of c-myc in the kidney relative to the brain samples and the final `relative_expression` of c-myc. In addition, an `error` term and a `lower` and `upper` intervals are provided.  

The previous analysis makes a few assumptions. One of which is a perfect amplification efficiency of the PCR reation. To assess the validity of this assumption, `pcr_assess` provides a method called `efficiency`. The input `data.frame` is the CT values of c-myc and GAPDH at different input amounts/dilutions.  

```r
## locate and read data
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- readr::read_csv(fl)

## make a vector of RNA amounts
amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

## calculate amplification efficiency
res <- pcr_assess(ct3,
                  amount = amount,
                  reference_gene = 'GAPDH',
                  method = 'efficiency')
res
```

In the case of using the Double Delta C_T, the assumption of the amplification efficiency is critical for the reliability of the model. In particulare, the *slope* and the *R^2* of the line between the log input amount and Delta C_T or differnce between the CT value of the target **c-myc** and GAPDH. Typically, The *slope* should be very small (less than 0.01). The `slope` here is appropriate, so the assumption holds true.  

### Other analysis methods `?pcr_analyze`  

* Delta CT method  
* Relative standard curve method  

### Testing statistical significance `?pcr_test`  

* Two group testing *t*-test and *wilcoxon* test  
* Linear regression testing  

## Documnetation  

```r
browseVignettes("pcr")
```
Alternatively, read the vingette can be found online, [here](http://rpubs.com/MahShaaban/pcr).

## Citation  

```r
citation("pcr")
```

## PeerJ Article 

For details about the methods and more examples, check out our [PeerJ article](https://peerj.com/articles/4473/?td=bl).

