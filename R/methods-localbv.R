#' @title Local Bivariate analysis with spatial autocorrelation
#' @description
#' This function is to explore the local bivariate relationship in the spatial space.
#' @rdname runLOCALBV-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param features1 the features name data object (only supporting character), see also \code{features2} parameter.
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
#' @param bv.method character one of the \code{'locallee'} and {'localmoran_bv'}, default is \code{'locallee'}.
#' @param bv.alternative a character string specifying the alternative hypothesis, default is \code{tow.sided}.
#' This only work when \code{bv.method = 'localmoran_bv'}.
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
#' @param lisa.method character one of the \code{'localG'} and {'localmoran'}, this is to perform the \code{LISA}
#' analysis using the result of bv.method, which can identify the spatial domain of the bivariate spatial
#' analysis result, default is \code{'localG'}.
#' @param lisa.alternative a character string specifying the alternative hypothesis, which works with
#' \code{lisa.method}, default is \code{greater}.
#' @param reduction.used character used as spatial coordinates to calculate the neighbours weights,
#' default is \code{NULL}, if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param permutation integer the permutation number to test, which only work with \code{bv.method='localmoran_bv'},
#' default is 100L.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis parallelly using
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel.
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)} to parallel it,
#' the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param action character, which control the type of return result, default is \code{get}, which will return
#' a \linkS4class{SimpleList}.
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
#' @return if \code{action = 'get'} (in default), the SimpleList object (like list object) will be return,
#' if \code{action = 'only'}, the data.frame will be return. if \code{action = 'add'}, the result of LISA is
#' stored in the \code{localResults} column of \code{int_colData} (internal column metadata). You can use
#' \code{localResults()} function of \code{SpatialFeatureExperiment} package to extract it.
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features, [`runGLOBALBV`] to
#' analysis the global bivariate spatial analysis, [`runLISA`] to identify the spatial domain of specified features.
#' @export
#' @author Shuangbin Xu
setGeneric('runLOCALBV',
  function(
    data,
    features1,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    bv.method = c("locallee", "localmoran_bv"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    reduction.used = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  )
  standardGeneric('runLOCALBV')
)


#' @rdname runLOCALBV-method
#' @aliases runLOCALBV,SingleCellExperiment
#' @export runLOCALBV
setMethod("runLOCALBV", "SingleCellExperiment", function(
    data,
    features1,
    features2 = NULL,
    assay.type = "logcounts",
    sample_id = 'all',
    bv.method = c("locallee", "localmoran"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    reduction.used = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  ){

  weight.method <- match.arg(weight.method)
  bv.method <- match.arg(bv.method)
  lisa.method <- match.arg(lisa.method)
  action <- match.arg(action)
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
                  n <- nrow(wm)
                  wi <- rowSums(wm)
                  wi2 <- rowSums(wm^2)
                  if (any(rowSums(wm) == 0)){
                      cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                  }
                  result <- .runLocalBv(xi, wm, features1, features2, n, permutation, bv.method, bv.alternative,
				       	random.seed, wi, wi2, lisa.method, lisa.alternative, BPPARAM)
                  return(result)
         })

  res <- .tidy_lisa_res(res)

  if (action == 'only'){
      res <- do.call("cbind", res)
      return(res)
  }

  if (action == 'add'){
      if (verbose){
          cli::cli_inform(c("The result is added to the input object, which can be extracted using",
                           paste0("`LISAResult()` with `type='", bv.method, ".SVP'`, and a specified `features`.")))
      }
      data <- .add_LISA_res(data, res, bv.method)
      return(data)
  }
  res <- S4Vectors::SimpleList(res)
  return(res)  
})


#' @rdname runLOCALBV-method
#' @aliases runLOCALBV,SVPExperiment
#' @export runLOCALBV
setMethod("runLOCALBV", "SVPExperiment",
  function(
    data,
    features1,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    bv.method = c("locallee", "localmoran_bv"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    reduction.used = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
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

       x2 <- .extract_gsvaExp_assay(data, gsvaexp, gsvaexp.assay.type)
       x2 <- x2[features2, , drop=FALSE]
       x <- rbind(x, x2)

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
                       n <- nrow(wm)
                       wi <- rowSums(wm)
                       wi2 <- rowSums(wm^2)
                       if (any(rowSums(wm) == 0)){
                           cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                       }
                       result <- .runLocalBv(xi, wm, features1, features2, n, permutation, bv.method, bv.alternative,
                                             random.seed, wi, wi2, lisa.method, lisa.alternative, BPPARAM)
                       return(result)
              })

       res <- .tidy_lisa_res(res)

       if (action == 'only'){
           res <- do.call("cbind", res)
           return(res)
       }

       if (action == 'add'){
           if (verbose){
               cli::cli_inform(c("The result is added to the input object, which can be extracted using",
                                paste0("`LISAResult()` with `type='", bv.method, ".SVP'`, and a specified `features`.")))
           }
           data <- .add_LISA_res(data, res, bv.method)
           return(data)
       }
       res <- S4Vectors::SimpleList(res)
    }else{
       res <- callNextMethod()
    }
    return(res)
})
