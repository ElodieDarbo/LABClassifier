
# LABClassifier

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of LABClassifier is to classify breast cancer samples into three subtypes: Luminal, Molecular Apocrine (MA) and Basal. To this end, it uses gene expression (RNA-seq, micro-array, illumina beadchip) of genes being determinant for the sensory/secretory split and Luminal/Apocrine split. We base our hypothesis on cell origin, specialisation and tumorogenesis.

## Installation

You can install the development version of LABClassifier from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ElodieDarbo/LABClassifier")
```

## Example

This is a basic example which shows you how to predict breast cancer subtypes from RNA-seq data (RSEM normalized):

``` r
library(LABClassifier)
## basic example code
data(TCGA.rsem)
LABclassifier(TCGA.rsem,plot=TRUE)
```

