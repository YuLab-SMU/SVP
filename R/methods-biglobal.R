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
#' @param random.seed numeric random seed number to repeatability, default is 1024.
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
    random.seed = 1024,
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
    random.seed = 1024,
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
                  res <- .internal.runBIGLOBAL(xi, wm, features1, features2, 
                                               permutation, alternative, 
                                               add.pvalue, random.seed)
                  return(res)
         })
  if (length(res) == 1){
      return(res[[1]])
  }
  names(res) <- sample_id
  res <- S4Vectors::SimpleList(res)
  return(res)
})

