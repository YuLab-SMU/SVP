context('fast_cor')

set.seed(123)
x <- matrix(rnorm(200), nrow=20)
x[x<0] <- 0
y <- Matrix::Matrix(x, sparse=T)

test_that('checking fast_cor',{
    r1 <- cor(t(x))
    res <- fast_cor(y)
    expect_equal(round(r1, 6), round(res$r, 6))
    r2 <- cor(t(x), method='spearman')
    res2 <- fast_cor(y, method = 'spearman')
    expect_equal(round(r2, 6), round(res2$r, 6))
})
