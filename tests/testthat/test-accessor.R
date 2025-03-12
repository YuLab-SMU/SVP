library(SingleCellExperiment) |> suppressPackageStartupMessages()
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
v <- log2(u + 1)
pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)

svpe <- SVPExperiment(assays=list(counts=u, logcounts=v), 
                     reducedDims=SimpleList(PCA=pca, tSNE=tsne))

test_that('using SVPExperiment to construct SVPExperiment object',{
  expect_true(inherits(svpe, "SVPExperiment"))
  expect_true(inherits(svpe, "SingleCellExperiment"))
})


sce1 <- SingleCellExperiment(matrix(rpois(1000, 5), ncol=ncol(svpe)))
rownames(sce1) <- paste0("GO:",seq(nrow(sce1)))
colnames(sce1) <- colnames(svpe)
sce2 <- SingleCellExperiment(matrix(rpois(1000, 5), ncol=ncol(svpe)))
rownames(sce2) <- paste0("KEGG:", seq(nrow(sce2)))
colnames(sce2) <- colnames(svpe)

gsvaExp(svpe, "GO") <- sce1
gsvaExp(svpe, "KEGG") <- sce2

test_that("using `gsvaExp<-` to set the value of gsvaExps in SVPExperiment",{
  elems <- int_colData(svpe)$gsvaExps 
  expect_true(length(elems)>0)

})


test_that("using gsvaExp to get the gsvaExp in SVPExperiment",{
  expect_true(identical(gsvaExp(svpe), sce1))
  expect_true(identical(gsvaExp(svpe, 2), sce2))
})


test_that("using gsvaExps to obtain all the gsvaExps in SVPExperiment",{
  elems <- gsvaExps(svpe)
  expect_equal(length(elems), 2)

})

test_that("using gsvaExpName to obtain the name of gsvaExp in SVPExperiment",{
  nm <- gsvaExpNames(svpe)
  expect_equal(c("GO", "KEGG"), nm)
})

test_that("`spatialCoords<-` to set the coordinate and spatialCoords to obtain the coordinate", {
  set.seed(123)
  da <- matrix(rnorm(200), ncol=2)
  spatialCoords(svpe) <- da
  expect_true(inherits(svpe, 'SVPExperiment'))
  
  dt <- spatialCoords(svpe)
  expect_true(inherits(dt, 'matrix') && ncol(dt)==2)
})

