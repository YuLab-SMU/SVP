set.seed(123)
coords <- matrix(abs(rnorm(100)), ncol=2)

weight <- SVP:::.obtain.weight(coords)

w2 <- deldir::deldir(coords) |>
     SVP:::.convert_to_distmt()

set.seed(123)
xx <- matrix(rnorm(500), ncol=50)
xx[xx<0] <- 0
xx <- Matrix::Matrix(xx, sparse=TRUE)

cal_lag <- function(x, y, w){
  res <- numeric(nrow(w))
  for (i in seq(ncol(w))){
    res[i] <- sum(x * w[i,]) * sum(y * w[i,])
  }
  return(res)
}

cal_lee <- function(x, y, W){
  n <- nrow(W)
  xx <- mean(x)
  yy <- mean(y)
  z <- x - xx
  zz <- sum(z^2)
  zy <- y - yy

  zzy <- sum(zy^2)

  llzy <- sum(cal_lag(z, zy, W))

  S2 <- sum(rowSums(W)^2)
  L <- (n/S2) * llzy/(sqrt(zz) * sqrt(zzy))
  return(L)
}

test_that("Lee calculation",{
  res1 <- SVP:::CalGlobalLeeParallel(xx, weight, f1=0, f2=1)
  res2 <- cal_lee(xx[1,], xx[2,], as.matrix(weight))
  expect_equal(as.numeric(res1$Lee), res2) 
}
)

cal_local_lee <- function(x, y, W){
  n <- nrow(W)
  xx <- mean(x)
  yy <- mean(y)

  z <- x - xx
  zy <- y - yy
  
  zz <- sum(z^2)
  zzy <- sum(zy^2)

  llzyi <- cal_lag(z, zy, W)

  Li <- (n * llzyi)/(sqrt(zz) * sqrt(zzy))
  
  return(Li)
}

test_that("Local lee calculation",{
  res1 <- SVP:::RunLocalLee(xx[1,], xx[2,], weight, n=nrow(weight))
  res2 <- cal_local_lee(xx[1,], xx[2,], as.matrix(weight))
  expect_equal(round(c(res1), 8), round(res2, 8))
})


data(hpda_spe_cell_dec)

test_that("runGLOBALBV test", {
  res1 <- runGLOBALBV(hpda_spe_cell_dec,
                      features1 = "Ductal APOL1 high-hypoxic",
                      features2 = c('Cancer clone A', "Cancer clone B"),
                      assay.type = 1
  )
  expect_true(inherits(res1, "list"))
  expect_equal(colnames(res1[[1]]), c('Cancer clone A', "Cancer clone B"))
  expect_equal(rownames(res1[[1]]), "Ductal APOL1 high-hypoxic")

})


test_that("runLOCALBV test", {
  res1 <- runLOCALBV(hpda_spe_cell_dec, 
                     features1 = "Ductal APOL1 high-hypoxic", 
                     features2 = c('Cancer clone A', "Cancer clone B"),
                     assay.type = 1
  ) 
  expect_true(inherits(res1, "SimpleList"))
})
