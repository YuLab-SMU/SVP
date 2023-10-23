<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVP: Efficient analysis of ‘gene set’ activity in spatial or single-cell data

SVP uses the KNN Graph and random walk with restrat based on the space
of MCA to calculate the activity of gene sets (such as kegg pathway,
signatures, go term, gene modules, …) efficiently. Then to detect the
spatially or single cell variable gene sets or (other features) with
avoiding the lost of low activity score of gene sets or other features
(such as the expression of genes), SVP uses the 2D Weighted Kernel
Density Estimation to process the score of gene sets (or expression of
genes) and uses Kullback–Leibler divergence to identify the signal
features.

# :writing_hand: Author

[Shuangbin Xu](https://github.com/xiangpin) and [Guangchuang
Yu](https://guangchuangyu.github.io)

School of Basic Medical Sciences, Southern Medical University

# :arrow_double_down: Installation

The development version from `github`:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("xiangpin/SVP")
```

# :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
