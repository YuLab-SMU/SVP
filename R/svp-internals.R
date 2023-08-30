#' @importFrom BiocParallel bplapply MulticoreParam
.cal_cor <- function(x, img, coords, beta = 2, cor.method = 'pearson', 
                     BPPARAM = NULL, threads = 5, p.adjust.method = 'fdr', ...){
    if (is.null(BPPARAM)) {
        BPPARAM <- MulticoreParam(workers = threads)
    }
    if (nrow(img) == 0){
        return(NULL)
    }else{
        coords <- coords * img[['scaleFactor']][1]
        img <- as.raster(img[['data']][[1]])
    }
    color.features <- .extract_color(img, beta = beta, coords)
    res <- bplapply(seq_len(nrow(x)), function(i){
        xi <- x[i,]
        tmp <- cor.test(as.numeric(xi), as.numeric(color.features), method = cor.method)
        res.i <- c(pvalue.cor.img = tmp$p.value, value.cor.img = as.numeric(tmp$estimate))
        return(res.i)
        }, 
        BPPARAM = BPPARAM
    )
    res <- do.call('rbind', res)
    p.adj.cor.img <- p.adjust(res[,1], method = p.adjust.method)
    res <- cbind(res, p.adj.cor.img)
    rownames(res) <- rownames(x)
    return(res)
}


