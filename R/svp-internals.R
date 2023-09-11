#' @importFrom spatialEco crossCorrelation
#' @importFrom BiocParallel bplapply MulticoreParam
.cal_cor <- function(x, img, coords, beta = 2, 
                     cor.method = 'pearson', 
                     BPPARAM = NULL, threads = 5, p.adjust.method = 'fdr', 
                     index.image = 1, ...){
    if (is.null(BPPARAM)) {
        BPPARAM <- MulticoreParam(workers = threads)
    }
    if (nrow(img) == 0){
        return(NULL)
    }else{
        if (is.numeric(index.image)){
            index.image <- unique(img$image_id)[index.image]
        }
        img <- img[img$image_id == index.image,]
        coords <- coords * img[['scaleFactor']]
        img <- as.raster(img[['data']][[1]])
    }
    coords.dist <- dist(coords)
    color.features <- .extract_color(img, beta = beta, coords)
    res <- bplapply(seq_len(nrow(x)), function(i){
        xi <- x[i,]
        #tmp <- cor.test(as.numeric(xi), as.numeric(color.features), method = cor.method)
        tmp <- crossCorrelation(as.numeric(xi), as.numeric(color.features), w = coords.dist)
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


