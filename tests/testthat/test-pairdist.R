context('pairdist')

set.seed(123)

a <- matrix(rnorm(60), nrow=10, ncol=6)

rownames(a) <- paste0('row', seq(nrow(a)))
colnames(a) <- paste0('col', seq(ncol(a)))

d <- pairDist(a, a)

test_that('checking pair distance ', {
  d2 <- dist(a) |> as.matrix()
  expect_equal(d, d2)
})
