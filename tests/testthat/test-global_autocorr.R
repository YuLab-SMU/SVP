set.seed(123)
coords <- matrix(abs(rnorm(100)), ncol=2)

weight <- SVP:::.obtain.weight(coords)

w2 <- deldir::deldir(coords) |> 
     SVP:::.convert_to_distmt()

test_that("weight matrix construction",{
  w3 <- SVP:::.norm_weight_mat(w2)
  expect_equal(as.matrix(weight)>0, as.matrix(w2)>0)
  expect_equal(weight, w3)

})

set.seed(123)
xx <- matrix(rnorm(500), ncol=50)
xx[xx<0] <- 0
xx <- Matrix::Matrix(xx, sparse=TRUE)

test_that("Moran's I calculation", {
  res1 <- SVP:::CalMoransiParallel(xx, wm=weight, permutation=1, lower_tail=0)
  res2 <- apply(as.matrix(xx), 1, 
                function(i){
                    ape::Moran.I(i, 
                                 weight = as.matrix(weight),
                                 alternative = 'greater'
                                 ) |> 
                      data.frame()})

  res2 <- do.call("rbind", res2)
  
  expect_equal(round(res1[, 1], 8), round(res2[, 1], 8))
  expect_equal(round(res1[, 2], 8), round(res2[, 2], 8))
  expect_equal(round(res1[, 3], 8), round(res2[, 3], 8))
  expect_equal(round(res1[, 5], 8), round(res2[, 4], 8))
})

cal_C <- function(x, W){
  n <- nrow(W)
  C <- (n - 1) * sum(W * (x - rep(x, each=n))^2) / (2 * sum(W) * sum((x - mean(x))^2))
  return(C)
}

test_that("Geary's C calculation", {
  res1 <- SVP:::CalGearyscParallel(xx, weight, permutation = 1, lower_tail = 1)
  res2 <- apply(as.matrix(xx), 1, cal_C, as.matrix(weight))
  expect_equal(round(res1[,1],8), round(res2, 8))
})

data(hpda_spe_cell_dec)

test_that("the result of runDetectSVG",{
  res <- hpda_spe_cell_dec |> runDetectSVG(assay.type = 1, verbose = FALSE)
  expect_true(length(svDfs(res))==1)
  expect_equal(rownames(svDf(res)), rownames(hpda_spe_cell_dec))
  expect_true(inherits(res, "SingleCellExperiment"))

  res <- res |> runDetectSVG(assay.type = 1, verbose = FALSE, method = 'geary')

  res <- res |> runDetectSVG(assay.type = 1, verbose = FALSE, method = 'getisord')
  expect_equal(c("sv.moransi", "sv.gearysc", "sv.getisord"), svDfNames(res))
})

test_that("the result of runKldSVG with kld", {
  res <- hpda_spe_cell_dec |> runKldSVG(assay.type = 1, verbose = FALSE)
  expect_true(inherits(svDf(res), "data.frame"))
})


