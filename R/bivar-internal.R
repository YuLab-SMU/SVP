.internal.runGLOBALBV <- function(
  x,
  weight,
  features1 = NULL,
  features2 = NULL,
  permutation = 100,
  alternative = c('greater', 'two.sided', 'less'),
  add.pvalue = FALSE,
  listn = NULL,
  across.gsvaexp = TRUE,
  random.seed = 1024
  ){
  allf <- rownames(x)
  alter <- switch(alternative, greater=3, `two.sided`=1, less = 2)
  if (is.null(features1) && is.null(features2)){
      f1 <- f2 <- seq(nrow(x))
  }else if (!is.null(features1) && is.null(features2)){
      f1 <- f2 <- .check_features(features1, allf, prefix='features1')
  }else if (is.null(features1) && !is.null(features2)){
      if (across.gsvaexp && length(listn) > 1){
          f1 <- .check_features(listn[[1]], allf, prefix='features2')
          f2 <- .check_features(unlist(listn[-1], use.names=FALSE), allf, prefix='features2')
      }else{
          f1 <- f2 <- .check_features(features2, allf, prefix="features2")
      }
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

.extract_gsvaExp_assay <- function(data, gsvaexp, gsvaexp.assay.type = NULL){
  if (is.null(gsvaexp.assay.type)){
      gsvaexp.assay.type <- 1
  }
  if (length(gsvaexp) > 1 && length(gsvaexp.assay.type)==1){
      gsvaexp.assay.type <- rep(gsvaexp.assay.type, length(gsvaexp))
  }
  res <- lapply(seq(length(gsvaexp)), function(x){
         assay(gsvaExp(data, gsvaexp[x]), gsvaexp.assay.type[x])
   })
  res <- do.call('rbind', res)
  return(res)
}

#' @importFrom utils combn
.runLocalBv <- function(
  x,
  weight,
  features1, 
  features2, 
  n, 
  listn = NULL,
  across.gsvaexp = TRUE,
  permutation = 200, 
  bv.method = c('locallee', 'localmoran_bv'),
  bv.alternative = c("two.sided", "greater", "less"),
  seed = 123,
  wi,
  wi2,
  lisa.method = c("localG", "localmoran"),
  lisa.alternative = c("greater", "two.sided", "less"),
  BPPARAM = SerialParam() 
){
  bv.method <- match.arg(bv.method)
  bv.alternative <- match.arg(bv.alternative)
  lisa.method <- match.arg(lisa.method)
  lisa.alternative <- match.arg(lisa.alternative)
  if (is.null(features1) && length(features2) > 1){
      if (length(listn)>1 && across.gsvaexp){
          allpair <- .generate_across_gsvaexp_pair(listn)
      }else{
          allpair <- combn(features2, 2) |> t()
      }
  }else if(is.null(features2) && length(features1) > 1){
      allpair <- combn(features1, 2) |> t()
  }else{
      allpair <- expand.grid(features1, features2) |> as.matrix()
  }
  if (bv.method == 'localmoran_bv'){
      res <- bplapply(seq(nrow(allpair)), function(i){
               .internal.runLocalMoranBv(x[allpair[i, 1],], x[allpair[i, 2],], weight, n, permutation, bv.alternative, 
                                      seed, wi, wi2, lisa.method, lisa.alternative)
               }, BPPARAM = BPPARAM)
  }else{
      res <- bplapply(seq(nrow(allpair)), function(i){
               .internal.runLocalLeeBv(x[allpair[i,1], ], x[allpair[i,2], ], weight, n, 
                                      wi, wi2, lisa.method, lisa.alternative)
               }, BPPARAM = BPPARAM)
  }
  names(res) <- paste(allpair[,1], allpair[,2], sep="_VS_")
  return(res)
}

.internal.runLocalMoranBv <- function(x, y, weight, n, permutation, bv.alternative, seed, 
                                      wi, wi2, lisa.method, lisa.alternative){
  res <- withr::with_seed(seed, RunLocalMoranBvPerm(x, y, weight, n, permutation)) |> 
         data.frame()
  prefix <- switch(bv.alternative, two.sided = "!=", greater = ">", less = "<")
  nm <- c("Ibvi", "E.Ibvi", "Var.Ibvi", "Z.Ibvi")
  pnm <- gettextf("Pr (z %s E(%s))", prefix, "Ibvi")
  if (bv.alternative == 'two.sided'){
      res[, 5] <- 2 * pnorm(abs(res[, 4]), lower.tail = FALSE)
  }else if (bv.alternative == 'greater'){
      res[, 5] <- pnorm(res[, 4], lower.tail = FALSE)
  }else{
      res[, 5] <- pnorm(res[, 4])
  }
  colnames(res) <- c(nm, pnm)
  lisa.res <- .internal.runLISA(res[, 1], weight, wi, wi2, n, lisa.method, lisa.alternative)
  res <- cbind(res, lisa.res)
  rownames(res) <- names(x)
  return(res)
}

.internal.runLocalLeeBv <- function(x, y, weight, n, wi, wi2, lisa.method, lisa.alternative){
  res <- RunLocalLee(x, y, weight, n) |> data.frame()
  colnames(res) <- "LocalLee"
  lisa.res <- .internal.runLISA(res[, 1], weight, wi, wi2, n, lisa.method, lisa.alternative)
  res <- cbind(res, lisa.res)
  rownames(res) <- names(x)
  return(res)
}

.generate_feature_listn <- function(x, f1, f2, ind){
  nm <- gsvaExpNames(x)
  if (is.numeric(ind)){
    ind <- nm[ind]
  }
  y <- lapply(ind, function(i) 
         intersect(f2, rownames(gsvaExp(x, i, withColData=FALSE, 
                                        withSpatialCoords=FALSE, withImgData=FALSE)))
  ) |> setNames(ind)
  if (is.null(f1)){
      return(y)
  }
  return(c(list(main = f1), y))
}

.tidy_globalbv_res <- function(x, y = NULL){
  lapply(x, function(i)as_tbl_df(i, y))
}

.generate_across_gsvaexp_pair <- function(x){
  pair1 <- names(x) |> combn(2) |> t()
  
  res <- lapply(seq(nrow(pair1)), function(i){
     x[pair1[i,]] |> expand.grid() |> as.matrix()
  })

  res <- do.call("rbind", res)

  return(res)

}

