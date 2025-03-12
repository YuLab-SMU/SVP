set.seed(123)
x <- matrix(rnorm(200), nrow=20)
rownames(x) <- paste0('r', seq(nrow(x)))
colnames(x) <- paste0('c', seq(ncol(x)))
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

    da <- as_tbl_df(res4)
    expect_true(inherits(da, "tbl_df"))
    expect_equal(dim(res4$r)[1] * dim(res4$r)[2], nrow(da))
})

data(hpda_spe_cell_dec)

test_that("runCORR test",{
    res1 <- runCORR(hpda_spe_cell_dec,
                    features1 = "Ductal APOL1 high-hypoxic",
                    features2 = c('Cancer clone A', "Cancer clone B"),
                    assay.type = 1,
                    action = 'only'
            )
    expect_true(inherits(res1, "tbl_df"))
})
