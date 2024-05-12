#' @title Bivariate analysis for spatial autocorrelation
#' @description
#' This function is to explore the bivariate relationship in the spatial space.
#' @rdname runBIGLOBAL-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param features1 the features name data object (only supporting character), default is NULL, 
#' see also \code{features2} parameter.
#' @param features2 character, if \code{features1} is not NULL, and \code{features2} is NULL, 
#' only the \code{features1} are analyzed, if \code{features1} is NULL, and \code{features2} is
#' is not NULL, the \code{features2} are analyzed, if \code{features2} is also NULL, all of features in the
#' \code{data} object will be analyzed. If \code{features2} and \code{features1} are not NULL, the bivariate
#' spatial autocorrelation analysis will be performed between the \code{features1} and \code{features2}.
#' default is \code{NULL}.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param sample_id character the sample(s) in the \linkS4class{SpatialExperiment} object whose cells/spots to use.
#' Can be \code{all} to compute metric for all samples; the metric is computed separately for each sample.
#' default is \code{"all"}.
#' @param method character now only the \code{'lee'}, default is \code{'lee'}.
#' @param weight object, which can be \code{nb}, \code{listw} or \code{Graph} object, default is NULL,
#' meaning the spatail neighbours weights will be calculated using the \code{weight.method}.
#' if the \code{data} contains multiple samples, and the \code{sample_id} is specified, it should be
#' provided as a list object with names (using \code{sample_id}).
#' @param weight.method character the method to build the spatial neighbours weights, default
#' is \code{knn} (k nearest neighbours). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is avaiable, such as \code{'tri2nb'}, \code{"knearneigh"},
#' \code{'dnearneigh'}, \code{"gabrielneigh"}, \code{"relativeneigh"}, which are from \code{spdep} package.
#' default is \code{knn}, if it is \code{"none"}, meaning the distance weight of each spot is used to
#' the weight.
#' @param reduction.used character used as spatial coordinates to calculate the neighbours weights, 
#' default is \code{UMAP}, if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param permutation integer the permutation number to test, default is 100L.
#' @param alternative a character string specifying the alternative hypothesis, which only work with 
#' \code{add.pvalue = TRUE}, default is \code{greater}.
#' @param add.pvalue logical whether calculate the pvalue, which is calculated with permutation test. So it might
#' be slow, default is \code{FALSE}, which the pvalue of result will be NULL.
#' @param verbose logical whether print the help information, default is TRUE.
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return SimpleList  
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features.
#' [`runLISA`] to explore the spatial hotspots.
#' @author Shuangbin Xu
#' @export
setGeneric('runBIGLOBAL',
  function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    method = c("lee"),
    weight = NULL,
    weight.method = c("knn", "tri2nb"),
    reduction.used = NULL,
    permutation = 100,
    alternative = 'greater',
    add.pvalue = FALSE,
    verbose = TRUE,
    ...
  )
  standardGeneric('runBIGLOBAL')
)

#' @rdname runBIGLOBAL-method
#' @aliases runBIGLOBAL,SingleCellExperiment
#' @export runBIGLOBAL
setMethod("runBIGLOBAL", "SingleCellExperiment", function(
    data, 
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    sample_id = 'all',
    method = c("lee"),
    weight = NULL, 
    weight.method = c("knn", "tri2nb"), 
    reduction.used = NULL,
    permutation = 100,
    alternative = 'greater',
    add.pvalue = FALSE,
    verbose = TRUE,
    ...
  ){
  
  weight.method <- match.arg(weight.method)
  method <- match.arg(method)
  sample_id <- .check_sample_id(data, sample_id)

  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)
  features <- union(features1, features2)

  if (length(features) >= 1){
      x <- x[features, ,drop=FALSE]
  }

  coords <- .check_coords(data, reduction.used, weight)

  res <- lapply(sample_id, function(sid){
                  if (sid == ".ALLCELL"){
                      ind <- seq(ncol(x))
                  }else{
                      ind <- colData(data)$sample_id == sid
                  }
                  coordsi <- if(!is.null(coords)){coords[ind, ,drop=FALSE]}else{NULL}
                  weighti <- if(inherits(weight, 'list')){weight[names(weight) == sid]}else{weight}
                  xi <- x[, ind, drop=FALSE]
                  wm <- .obtain.weight(coordsi, weight = weighti, weight.method = weight.method, ...)
                  if (any(rowSums(wm) == 0)){
                      cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                  }
                  res <- .internal.runBIGLOBAL(xi, wm, features1, features2, permutation, alternative, add.pvalue)
                  return(res)
         })
  if (length(res) == 1){
      return(res[[1]])
  }
  res <- S4Vectors::SimpleList(res)
  return(res)
})


.internal.runBIGLOBAL <- function(
  x, 
  weight,
  features1 = NULL, 
  features2 = NULL, 
  permutation = 100, 
  alternative = c('greater', 'two.sided', 'less'), 
  add.pvalue = FALSE
  ){
  allf <- rownames(x)
  alter <- switch(alternative, greater=2, `two.sided`=1, less = 3)
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
  
  L <- CalGlobalLeeParallel(x, weight, f1-1L, f2-1L, permutation, alter, add.pvalue)
  rownames(L) <- allf[f1]
  colnames(L) <- allf[f2]
  if (add.pvalue){
      pv <- CalGlobalLeeParallel(x, weight, f1-1L, f2-1L, permutation, alter, add.pvalue)
      rownames(pv) <- allf[f1]
      rownames(pv) <- allf[f2]
  }else{
      pv <- NULL
  }
  return(list(Lee=L, pvalue=pv))
}


.check_features <- function(x, y, prefix){
  x <- unique(x)
  f1 <- match(x, y)
  f1 <- f1[!is.na(f1)]
  if (length(f1) < 1){
      cli::cli_abort(paste0("The `", prefix[1],"` is/are not present in the row names."))
  }
  return(f1)
}


.check_coords <- function(data, reduction.used, weight = NULL){
  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(reduction.used %in% reducedDimNames(data))
  coords <- NULL
  if((flag1 || flag2) && is.null(weight)){
      if (flag2){
          coords <- reducedDim(data, reduction.used)
          coords <- coords[,c(1, 2)]
      }
      if (flag1 ){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
  }else if (all(!flag1, !flag2) && is.null(weight)){
      cli::cli_abort(c("The {.cls {class(data)}} should have 'spatialCoords' or the reduction result of 'UMAP' or 'TSNE'.
                     Or the `weight` should be provided."))
  }

  return(coords)
}
