# Inference-for-Dispersion-and-Curvature

This repository contains the implementation for the paper "Inference for Dispersion and Curvature of Random Objects" in R.

### Supporting software requirements

R version 4.3.1.

### Libraries and dependencies used by the code

R package to run simulations, real data analysis and to reproduce figures:

* igraph version 1.5.0
* ellipse version 0.4.5
* CovTools version 0.5.4
* spherepc version 0.1.7
* abind version 1.4-5
* fdapace version 0.5.9
* colorspace version 2.1-0
* plotrix version 3.8-2
* expm version 0.999-7
* vegan version 2.6-4
* readxl version 1.4.3
* fdadensity version 0.1.2
* ggplot2 version 3.5.0
* fields version 14.1
* scatterplot3d version 0.3-44
* expm version 0.999-7

All packages are available through CRAN (https://cran.r-project.org/) and can be installed automatically by running install.packages(PACKAGE_NAME) in an R session.

### Folder Structure

* `./src/`  code for all functions used in the paper.
* `./simulations/`  code to reproduce simulations in Section 5 and Supplment Section S.1.
* `./Applications/` code to reproduce data analysis in Section 6 and Supplment Section S.2.
