set.seed(123)
coords <- matrix(abs(rnorm(100)), ncol=2)

weight <- SVP:::.obtain.weight(coords)

set.seed(123)
xx <- matrix(rnorm(500), ncol=50)
xx[xx<0] <- 0
xx <- Matrix::Matrix(xx, sparse=TRUE)


cal_localmoran <- function(x, W){
  n <- nrow(W)
  local_moran <- numeric(n)  
  x_bar <- mean(x)  
  s2 <- sum((x - x_bar)^2) / n  
  
  for (i in 1:n) {  
    z_i <- x[i] - x_bar  
  
    w_ij_z_j <- sum(W[i, ] * (x - x_bar))  
  
    local_moran[i] <- (z_i / s2) * w_ij_z_j  
  } 
  return(local_moran) 
}


test_that("Local Moran calculation", {
  res1 <- SVP:::CalLocalMoranParallel(xx[1,,drop=FALSE], weight)
  res2 <- cal_localmoran(xx[1,], weight)
  expect_equal(round(res1[[1]][,1], 8), round(res2, 8))
})

data(hpda_spe_cell_dec)

test_that("runLISA test and LISAsce test",{
  res <- hpda_spe_cell_dec |> runLISA(assay.type=1, features=rownames(hpda_spe_cell_dec))
  
  res2 <- LISAsce(hpda_spe_cell_dec, res)

  expect_true(inherits(res, 'SimpleList'))
  expect_equal(rownames(hpda_spe_cell_dec), names(res))
  expect_true(inherits(res2, "SVPExperiment"))
  
})

test_that("runLISA test with action = 'add' and  LISAResult to obtain the result",{
  res <- hpda_spe_cell_dec |> 
      runLISA(assay.type = 1, features = rownames(hpda_spe_cell_dec), action = 'add')
  expect_true(inherits(res, "SingleCellExperiment"))

  res2 <- LISAResult(res, type='localG.SVP', features=rownames(hpda_spe_cell_dec))
  expect_true(inherits(res2, 'SimpleList'))
  expect_equal(rownames(hpda_spe_cell_dec), names(res2))
})
