<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVP: Efficient analysis of ‘gene set’ activity in spatial or single-cell data

## :newspaper: Description

SVP uses the distance between cells and cells, features and features,
cells and features in the space of MCA to build nearest neighbor graph,
then uses random walk with restart algorithm to calculate the activity
score of gene sets (such as kegg pathway, signatures, go term, gene
modules, transcription factor, …), which is then further weighted using
the hypergeometric test results from the original expression matrix. In
addition, to detect the spatially or single cell variable gene sets or
(other features) accurately, SVP uses the 2d weighted kernel density
estimation to process the score of gene sets (or expression of genes)
and uses Kullback–Leibler divergence to identify the spatial variable
features based on permutation test. SVP also provides Geary’s and
Morans’I that measure spatial autocorrelation to identify the spatial
variable features efficiently based on Rcpp and RcppParallel. SVP is
developed based on SingleCellExperiment class, which can be
interoperable with the existing computing ecosystem.

## :writing_hand: Author

[Shuangbin Xu](https://github.com/xiangpin) and [Guangchuang
Yu](https://guangchuangyu.github.io)

School of Basic Medical Sciences, Southern Medical University

## :arrow_double_down: Installation

The development version from `github`:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("xiangpin/SVP")
```

To enhance performance, it is **strongly recommended** to connect your R
BLAS library with the
[OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) library for matrix
calculations. This can be accomplished using the
[ropenblas](https://prdm0.github.io/ropenblas/) package. Or you can
install [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) and link the
library to R library by
`ln -s your_openblas_installed_path_libopenblas.so your_R_install_path_libRblas.so`
manually.

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
