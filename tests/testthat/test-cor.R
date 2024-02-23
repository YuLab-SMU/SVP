context('fast_cor')

set.seed(123)
x <- matrix(rnorm(200), nrow=20)
x[x<0] <- 0
x1 <- x[seq(8),]
x2 <- x[seq(9,20),]
y <- Matrix::Matrix(x, sparse = TRUE)
y1 <- Matrix::Matrix(x1, sparse = TRUE)
y2 <- Matrix::Matrix(x2, sparse = TRUE)

test_that('checking fast_cor',{
    r1 <- cor(t(x))
    res <- fast_cor(y)
    expect_equal(round(r1, 6), round(res$r, 6))
    r2 <- cor(t(x), method='spearman')
    res2 <- fast_cor(y, method = 'spearman')
    expect_equal(round(r2, 6), round(res2$r, 6))
    r3 <- cor(t(x1), t(x2))
    res3 <- fast_cor(y1, y2)
    expect_equal(round(r3, 6), round(res3$r, 6))
    r4 <- cor(t(x1), t(x2), method = 'spearman')
    res4 <- fast_cor(y1, y2, method = 'spearman')
    expect_equal(round(r4, 6), round(res4$r, 6))

})
