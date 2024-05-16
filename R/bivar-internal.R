.internal.runGLOBALBI <- function(
  x,
  weight,
  features1 = NULL,
  features2 = NULL,
  permutation = 100,
  alternative = c('greater', 'two.sided', 'less'),
  add.pvalue = FALSE,
  random.seed = 1024
  ){
  allf <- rownames(x)
  alter <- switch(alternative, greater=3, `two.sided`=1, less = 2)
  if (is.null(features1) && is.null(features2)){
      f1 <- f2 <- seq(nrow(x))
  }else if (!is.null(features1) && is.null(features2)){
      f1 <- f2 <- .check_features(features1, allf, prefix='features1')
  }else if (is.null(features1) && !is.null(features2)){
      f1 <- f2 <- .check_features(features2, allf, prefix="features2")
  }else if (!is.null(features1) && !is.null(features2)){
      f1 <- .check_features(features1, allf, prefix='features1')
      f2 <- .check_features(features2, allf, prefix='features2')
  }

  res <- withr::with_seed(random.seed, 
         CalGlobalLeeParallel(x, weight, f1-1L, f2-1L, 
                              permutation, alter, add.pvalue))
  rownames(res$Lee) <- allf[f1]
  colnames(res$Lee) <- allf[f2]
  if (add.pvalue){
      rownames(res$pvalue) <- allf[f1]
      colnames(res$pvalue) <- allf[f2]
  }else{
      res$pvalue <- NULL
      res <- c(res, list(pvalue=NULL))
  }
  return(res)
}

