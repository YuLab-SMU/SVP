.internal.runBIGLOBAL <- function(
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

  L <- CalGlobalLeeParallel(x, weight, f1-1L, f2-1L, permutation, alter, FALSE)
  rownames(L) <- allf[f1]
  colnames(L) <- allf[f2]
  if (add.pvalue){
      pv <- withr::with_seed(random.seed,
               CalGlobalLeeParallel(x, weight, f1-1L, f2-1L, permutation, alter, TRUE)
      )
      rownames(pv) <- allf[f1]
      colnames(pv) <- allf[f2]
  }else{
      pv <- NULL
  }
  return(list(Lee=L, pvalue=pv))
}

