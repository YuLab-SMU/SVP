#' @title Bivariate analysis for spatial autocorrelation
#' @description
#' This function is to explore the bivariate relationship in the spatial space.
#' @rdname runGLOBALBI-method
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
#' is \code{voronoi} (Voronoi tessellation). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is avaiable, such as \code{"knearneigh"},
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
#' @param gsvaexp character the one character from the name of `gsvaExpNames(data)`, default is NULL. If \code{data}
#' is \linkS4class{SVPExperiment}, and the parameter is specified simultaneously. the `features` (Usually genes)
#' from the displayed class, and `gsvaexp.features` from name in `rownames(gsvaExp(data, gsvaexp))` will be
#' performed the analysis.
#' @param gsvaexp.assay.type character the assay name in the `assays(gsvaExp(data, gsvaexp))`, default is NULL, 
#' which works with \code{gsvaexp} parameter.
#' @param gsvaexp.features character the name from the `rownames(gsvaExp(data, gsvaexp))`. If \code{gsvaexp} is
#' specified and \code{data} is \linkS4class{SVPExperiment}, it should be provided. Default is NULL.
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return SimpleList  
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features.
#' [`runLISA`] to explore the spatial hotspots.
#' @author Shuangbin Xu
#' @export
setGeneric('runGLOBALBI',
  function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    method = c("lee"),
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    reduction.used = NULL,
    permutation = 100,
    alternative = 'greater',
    add.pvalue = FALSE,
    random.seed = 1024,
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,    
    ...
  )
  standardGeneric('runGLOBALBI')
)

#' @rdname runGLOBALBI-method
#' @aliases runGLOBALBI,SingleCellExperiment
#' @export runGLOBALBI
setMethod("runGLOBALBI", "SingleCellExperiment", function(
    data, 
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    sample_id = 'all',
    method = c("lee"),
    weight = NULL, 
    weight.method = c("voronoi", "knn", "none"), 
    reduction.used = NULL,
    permutation = 100,
    alternative = 'greater',
    add.pvalue = FALSE,
    random.seed = 1024,
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,    
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
                  res <- .internal.runGLOBALBI(xi, wm, features1, features2, 
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


#' @rdname runGLOBALBI-method
#' @aliases runGLOBALBI,SVPExperiment
#' @export runGLOBALBI
setMethod("runGLOBALBI", "SVPExperiment", function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    sample_id = 'all',
    method = c("lee"),
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    reduction.used = NULL,
    permutation = 100,
    alternative = 'greater',
    add.pvalue = FALSE,
    random.seed = 1024,
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  ){

    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform(c("The {.var gsvaexp} is specified, the specified {.var gsvaExp} will be",
                           "used to perform the analysis with `features` displayed from `rownames(data)`."))
       }
       weight.method <- match.arg(weight.method)
       method <- match.arg(method)
       sample_id <- .check_sample_id(data, sample_id)

       if (is.null(assay.type)){
           assay.type <- 1
       }

       x <- assay(data, assay.type)
       if (is.null(features1) && is.null(features2)){
           cli::cli_abort(c("The {.var gsvaexp} is specified, and the `data` is {.cls {class(data)}}.",
                            "The `features1` and `features2` should not be `NULL` simultaneously."))
       }
       if (!is.null(features1) && !is.null(features2)){
           cli::cli_inform(c("The `data` is {.cls {class(data)}} class, the `features1` and `features2` which are from ",
                            "the displayed object, will be merged to perform analysis with `gsvaexp.features`."))
       }
       if (is.null(gsvaexp.features) && !is.null(gsvaexp)){
           cli::cli_abort(c(paste0("The `data` is {.cls {class(data)}} class and `gsvaexp='", gsvaexp,"'`."),
                          " The `gsvaexp.features` should be provided, you can use ",
                          "`rownames(gsvaExp(data)) |> head()` to view them."))
       }

       features1 <- union(features1, features2)
       features2 <- gsvaexp.features

       if (length(features1) >= 1){
           x <- x[features1, , drop=FALSE]
       }

       data2 <- gsvaExp(data, gsvaexp)
       if (is.null(gsvaexp.assay.type)){
           gsvaexp.assay.type <- 1
       }
       x2 <- assay(data2, gsvaexp.assay.type)
       x2 <- x2[features2, , drop=FALSE]
       x <- Matrix::rbind2(x, x2)

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
                       res <- .internal.runGLOBALBI(xi, wm, features1, features2,
                                                    permutation, alternative,
                                                    add.pvalue, random.seed)
                       return(res)
              })
       names(res) <- sample_id
       res <- S4Vectors::SimpleList(res)
    }else{
       res <- callNextMethod()
    }
    return(res)
})
