<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVP: Robust analysis of ‘gene set’ activity in spatial or single-cell data

## :newspaper: Description

SVP uses the distance between cells and cells, features and features,
cells and features in the space of MCA to build nearest neighbor graph,
then uses random walk with restart algorithm to calculate the activity
score of gene sets (such as cell marker genes, kegg pathway, go
ontology, gene modules, transcription factor or miRNA target sets,
reactome pathway, …), which is then further weighted using the
hypergeometric test results from the original expression matrix. To
detect the spatially or single cell variable gene sets or (other
features) and the spatial colocalization between the features
accurately, SVP provides some global and local spatial autocorrelation
method to identify the spatial variable features. SVP is developed based
on SingleCellExperiment class, which can be interoperable with the
existing computing ecosystem.

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
library to R BLAS library by
`ln -s your_openblas_installed_path_libopenblas.so your_R_install_path_libRblas.so`
manually.

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
