# scFeatureFilter
By Angeles Arzalluz-Luque, Guillaume Devailly, Anna Mantsoki & Anagha Joshi.

![scFeatureFilter outputs](figure1_small.png)


scFeatureFilter is an R package that assists the identification and removal of undesired technical variability in single-cell RNAseq datasets. Our tool provides functions to explore the technical and biological variability components in the expression data, aiming to help set a threshold for the selection of features (genes, transcripts...) that are lowly affected by noise (i.e. feature filtering). Ultimately, this will help obtain more reliable biological conclusions from single-cell RNA-seq data.

# Install
At the moment, only the development version of scFeatureFilter can be installed. Use the `install_github` function in the `devtools` package (requires R version ≥ 3.4.2) as indicated below:
```R
devtools::install_github("gdevailly/scFeatureFilter")
```

# Getting started
Load the package:
```R
library(scFeatureFilter)
```

Start by reading the package vignette:
```R
vignette("Introduction", package = "scFeatureFilter")
```
