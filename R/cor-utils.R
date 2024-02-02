#' Calculation of correlations and associated p-values
#' @param x sparse Matrix
#' @param y sparse Matrix which has the same column length of \code{x},
#' default is NULL.
#' @param combine logical whether combine the correlation of \code{x} and \code{y}
#' if \code{y} is provided, default is FALSE.
#' @param method a character string indicating which correlation coefficient,
#' One of \code{"pearson"} (default), \code{"spearman"} and \code{"bicorr"}.
#' @param alternative indicates the alternative hypothesis and must be one of
#' the initial letter \code{"two.sided"}, \code{"less"} and \code{"greater"}.  
#' \code{"greater"} corresponds to positive association, \code{"less"} to negative 
#' association, default is \code{"two.sided"}.
#' @return a list containing the matrix of correlation and matrix of pvalue.
#' @importFrom stats pt
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol=10) 
#' x <- Matrix::Matrix(x, sparse = TRUE)
#' x1 <- x[seq(10),]
#' x2 <- x[seq(11, 50),]
#' res <- fast_cor(x = x1, y = x2, combine = FALSE)
#' res$r |> dim()
#' res2 <- fast_cor(x = x1, y = x2, combine = TRUE)
#' res2$r |> dim()
fast_cor <- function(x, 
                     y = NULL, 
                     combine = FALSE,
                     method = c('pearson', 'spearman', 'bicorr'), 
                     alternative = c("two.sided", "less", "greater")
                     ){

    method <- match.arg(method)
    ia = match.arg(alternative)
    ny <- 0
    indx <- !is.na(x)
    np <- NULL 
    if (!is.null(y)){
        indy <- !is.na(y)
        x <- Matrix::rbind2(x, y)
        if (!combine){
            indy <- !is.na(y)
            ny <- nrow(y)
            np <- indx %*% Matrix::t(indy) |> as.matrix()
        }
    }
        
    if (is.null(np)){
        indx <- !is.na(x)
        np <- indx %*% Matrix::t(indx) |> as.matrix()
    }
    nx <- nrow(x)
    nx <- nx - ny
    if (method == 'spearman'){
        x <- DelayedMatrixStats::rowRanks(x, ties.method = 'average', useNames=TRUE) 
        x <- Matrix::Matrix(x, sparse=TRUE)
    }
    if (method == 'bicorr'){
        mc <- CalParallelBiCor(x)
    }else{
        mc <- CalParallelCor(x)
    }

    rownames(mc) <- colnames(mc) <- rownames(x)

    if (!is.null(y) && nx != nrow(x)){
        mc <- mc[seq(nx), seq(nx+1, nx+ny)]
    }
    lower.tail <- FALSE
    mc2 <- mc
    if (ia == "two.sided") {
        mc2 <- abs(mc)
    }else if (ia == 'less'){
        lower.tail <- TRUE
    }
    t.val = sqrt(np - 2) * mc2/sqrt(1 - mc2^2)
    p = pt(t.val, np - 2, lower.tail = lower.tail)
    if (ia == 'two.sided'){
        p <- 2 * p
    }

    return(list(r=mc, pval=p))
}


