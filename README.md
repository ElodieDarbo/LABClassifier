
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
# Import espression matrices from Zenodo doi:10.5281/zenodo.10935179
# dataset is a character vector containing one or more dataset identifiers
# from "TCGA", "METABRIC", "ICGC", "EORTC_x3p" or "all" if all four datasets are wanted.
dataset <- "TCGA"
data <- import.data.zenodo(dataset = dataset)
# If the data were already dowloaded, they are stored in "./ext_data" and can
# be accessed using read.table("./ext_data/paste(dataset,"matrix.txt",sep="_"),header=T,row.names=1)
# Test on TCGA data
predictions <- LABclassifier(data[[dataset]],plot=TRUE)
```

## References

Deena M.A. Gendoo, Natchar Ratanasirigulchai, Markus S. Schroeder, Laia Pare, Joel S Parker, Aleix Prat, Benjamin Haibe-Kains (2023). genefu: Computation of Gene Expression-Based Signatures in Breast Cancer. [doi:10.18129/B9.bioc.genefu](https://doi.org/10.18129/B9.bioc.genefu), [R package version 2.34.0](https://bioconductor.org/packages/genefu).

Hänzelmann, S., Castelo, R. and Guinney, A. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics, 14:7, 2013.

Scrucca L, Fraley C, Murphy TB, Raftery AE (2023). Model-Based Clustering, Classification, and Density Estimation Using mclust in R. Chapman and Hall/CRC. ISBN 978-1032234953, [doi:10.1201/9781003277965](https://doi.org/10.1201/9781003277965), [book](https://mclust-org.github.io/book/).

Kolde R (2019). pheatmap: Pretty Heatmaps. [R package version 1.0.12](https://CRAN.R-project.org/package=pheatmap).

Blondel E (2023). zen4R: Interface to 'Zenodo' REST API. [R package version 0.9](https://CRAN.R-project.org/package=zen4R).

H. Wickham (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2024). viridis(Lite) - Colorblind-Friendly Color Maps for R. viridis package version 0.6.5.
