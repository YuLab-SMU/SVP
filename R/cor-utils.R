#' Calculation of correlations and associated p-values
#' @param x sparse Matrix which rows are the features and columns are the samples.
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
#' @param add.pvalue logical whether calculate the pvalue of correlation using t 
#' test, default is FALSE.
#' @return a list containing the matrix of correlation and matrix of pvalue (if 
#' \code{add.pvalue} is FALSE (default), the matrix of pvalue will be NULL).
#' @importFrom stats pt
#' @export
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol=10)
#' rownames(x) <- paste0('row', seq(nrow(x)))
#' colnames(x) <- paste0('col', seq(ncol(x)))
#' x <- Matrix::Matrix(x, sparse = TRUE)
#' x1 <- x[seq(10),]
#' x2 <- x[seq(11, 50),]
#' res <- fast_cor(x = x1, y = x2, combine = FALSE)
#' res$r |> dim()
#' res2 <- fast_cor(x = x1, y = x2, combine = TRUE)
#' res2$r |> dim()
fast_cor <- function(
    x, 
    y = NULL, 
    combine = FALSE,
    method = c('pearson', 'spearman', 'bicorr'), 
    alternative = c("two.sided", "less", "greater"),
    add.pvalue = FALSE
){
    if (!inherits(x, 'dgCMatrix')){
        x <- suppressWarnings(Matrix::Matrix(as.matrix(x), sparse=TRUE))
    }
    method <- match.arg(method)
    ia <- match.arg(alternative)
    #indx <- !is.na(x)
    np <- NULL
    if (!is.null(y)){
        if (ncol(x) != ncol(y)){
            cli::cli_abort("The column number of {.var x} and {.var y} should be equal when {.var y} is provided!")
        }
        if (!inherits(y, 'dgCMatrix')){
            y <- suppressWarnings(Matrix::Matrix(as.matrix(y), sparse=TRUE))
        }
        if (combine){
            x <- Matrix::rbind2(x, y)
        }else{
            if (add.pvalue){
                #indy <- !is.na(y)
                #np <- indx %*% Matrix::t(indy) |> as.matrix()
                np <- .mat_mult(x, y)
            }
        }
    }
        
    if (is.null(np) && add.pvalue){
        #indx <- !is.na(x)
        #np <- indx %*% Matrix::t(indx) |> as.matrix()
        np <- .mat_mult(x)
    }

    if (method == 'spearman'){
        x <- DelayedMatrixStats::rowRanks(x, ties.method = 'average', useNames=TRUE) 
        x <- Matrix::Matrix(x, sparse=TRUE)
        if (!is.null(y) && !combine){
            y <- DelayedMatrixStats::rowRanks(y, ties.method = 'average', useNames=TRUE)
            y <- Matrix::Matrix(y, sparse=TRUE)
        }
    }

    if (method == 'bicorr'){
        if (!is.null(y) && !combine){
            mc <- CalParallelBiCorTwoMatrix(x, y)
            rownames(mc) <- rownames(x)
            colnames(mc) <- rownames(y) 
        }else{
            mc <- CalParallelBiCor(x)
            rownames(mc) <- colnames(mc) <- rownames(x)
        }
    }else{
        if (!is.null(y) && !combine){
            mc <- corCpp(Matrix::t(x), Matrix::t(y))
            rownames(mc) <- rownames(x)
            colnames(mc) <- rownames(y)
        }else{
            mc <- CalParallelCor(x)
            rownames(mc) <- colnames(mc) <- rownames(x)
        }
    }

    if (add.pvalue){
        p <- .cal_cor_p(mc, ia, np)
    }else{
        p <- NULL
    }
    
    return(list(r=mc, pval=p))
}

.cal_cor_p <- function(mc, ia, np){
    lower.tail <- FALSE
    if (ia == "two.sided"){
        mc <- abs(mc)
    }else if (ia == 'less'){
        lower.tail <- TRUE
    }
    t.val = sqrt(np - 2) * mc/sqrt(1 - mc^2)
    p = pt(t.val, np - 2, lower.tail = lower.tail)
    if (ia == 'two.sided'){
        p <- 2 * p
    }

    return(p)
}

.mat_mult <- function(x, y = NULL){
    x <- .convert_flag(x)
    if (!is.null(y)){
        y <- .convert_flag(y)
        res <- MatMultCpp(x, t(y))
    }else{
        res <- MatMultCpp(x, t(x))
    }
    return(res)
}

.convert_flag <- function(x){
    x <- as.matrix(!is.na(x))
    x[x] <- 1
    x[!x] <- 0
    return(x)
}


.internal.runCORR <- function(
    x, 
    features1, 
    features2, 
    method = "spearman",
    alternative = 'two.sided',
    add.pvalue = FALSE,
    listn = NULL,
    across.gsvaexp = TRUE
){
  allf <- rownames(x)
  y <- NULL
  if (is.null(features1) && !is.null(features2)){
      if (across.gsvaexp && length(listn) > 1){
          f1 <- .check_features(listn[[1]], allf, prefix='features2')
          f2 <- .check_features(unlist(listn[-1], use.names=FALSE), allf, prefix='features2')
          y <- x[f2, ,drop=FALSE]
          x <- x[f1, ,drop=FALSE]
      }else{
          f1 <- .check_features(features2, allf, prefix="features2")
          x <- x[f1,,drop=FALSE]
      }
  }else if (!is.null(features1) && !is.null(features2)){
      f1 <- .check_features(features1, allf, prefix='features1')
      f2 <- .check_features(features2, allf, prefix='features2')
      y <- x[f2,,drop=FALSE]
      x <- x[f1,,drop=FALSE]
  }
  res <- fast_cor(x = x, 
                  y = y, 
                  method = method, 
                  alternative = alternative, 
                  add.pvalue = add.pvalue
  )
  return(res)
}
