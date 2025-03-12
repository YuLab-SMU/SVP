#' @title Local Bivariate analysis with spatial autocorrelation
#' @description
#' This function is to explore the local bivariate relationship in the spatial space. Like
#' \code{runGLOBALBV}, It efficiently reflects the extent to which bivariate associations 
#' are spatially grouped in local. Put differently, it can be utilized to quantify the 
#' bivariate spatial dependency in local. See also the references.  
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
#' @param bv.method character one of the \code{'locallee'} and \code{'localmoran_bv'}, default is \code{'locallee'}.
#' @param bv.alternative a character string specifying the alternative hypothesis, default is \code{tow.sided}.
#' This only work when \code{bv.method = 'localmoran_bv'}.
#' @param weight object, which can be \code{nb}, \code{listw} or \code{Graph} object, default is NULL,
#' meaning the spatial neighbours weights will be calculated using the \code{weight.method}.
#' if the \code{data} contains multiple samples, and the \code{sample_id} is specified, it should be
#' provided as a list object with names (using \code{sample_id}).
#' @param weight.method character the method to build the spatial neighbours weights, default
#' is \code{voronoi} (Voronoi tessellation). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is available, such as \code{"knearneigh"},
#' \code{'dnearneigh'}, \code{"gabrielneigh"}, \code{"relativeneigh"}, which are from \code{spdep} package.
#' default is \code{knn}, if it is \code{"none"}, meaning the distance weight of each spot is used to
#' the weight.
#' @param lisa.method character one of the \code{'localG'} and \code{'localmoran'}, this is to perform the \code{LISA}
#' analysis using the result of bv.method, which can identify the spatial domain of the bivariate spatial
#' analysis result, default is \code{'localG'}.
#' @param lisa.alternative a character string specifying the alternative hypothesis, which works with
#' \code{lisa.method}, default is \code{greater}.
#' @param lisa.flag.method a character string specifying the method to calculate the threshold for the cluster
#' type, default is \code{"mean"}. Other option is \code{"median"}. 
#' @param reduction.used character used as spatial coordinates to calculate the neighbours weights,
#' default is \code{NULL}, the result of reduction can be specified, such as \code{UMAP}, \code{TSNE}, \code{PCA}.
#' If it is specified, the weight neighbours matrix will be calculated using the result of specified reduction.
#' @param group.by character a specified category column names (for example the cluster column name) of
#' \code{colData(data)}. Or a vector of length equal to ‘ncol(x)’, specifying the group to which each cell
#' is assigned. If it was specified, the adjacency weighted matrix will be built based on the principle that
#' spots or cells in the same category are adjacent, default is NULL.
#' @param permutation integer the permutation number to test, which only work with \code{bv.method='localmoran_bv'},
#' default is 100L.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis in parallel using
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel.
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=TRUE)} to parallel it,
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
#' @param across.gsvaexp logical whether only calculate the relationship of features between the multiple `gsvaExps`
#' not the internal features of gsvaExp. For example, \code{'a'} and \code{'b'} features are from the \code{'AB'} `gsvaExp`,
#' \code{'c'} and \code{'d'} features are from the \code{'CD'} `gsvaExp`. When \code{across.gsvaexp=TRUE} and
#' \code{gsvaexp.features = c('a', 'b', 'c', 'd')} and \code{gsvaexp = c('AB', 'CD')}, Only the relationship of
#' \code{a} and \code{c}, \code{a} and \code{d}, \code{b} and \code{c}, and \code{b} and \code{d} will be calculated.
#' default is TRUE.
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return if \code{action = 'get'} (in default), the SimpleList object (like list object) will be return,
#' if \code{action = 'only'}, the data.frame will be return. if \code{action = 'add'}, the result of LISA is
#' stored in the \code{localResults} column of \code{int_colData} (internal column metadata). You can use
#' \code{localResults()} function of \code{SpatialFeatureExperiment} package to extract it.
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features, [`runGLOBALBV`] to
#' analysis the global bivariate spatial analysis, [`runLISA`] to identify the spatial domain of specified features.
#' @export
#' @references
#' Lee, SI. Developing a bivariate spatial association measure: An integration of Pearson's r and Moran's I . 
#' J Geograph Syst 3, 369–385 (2001). https://doi.org/10.1007/s101090100064
#' @author Shuangbin Xu
#' @examples
#' data(hpda_spe_cell_dec)
#' res1 <- hpda_spe_cell_dec |> runLOCALBV(
#'           features1 = 'Cancer clone A', 
#'           features2 = 'Cancer clone B', 
#'           assay.type='affi.score'
#'         )
#' res1
#' res1[['Cancer clone A_VS_Cancer clone B']] |> head()
#' # add the LocalLee and Gi of LOCALBV result to input object.
#' hpda_spe_cell_dec <- LISAsce(hpda_spe_cell_dec, res1, 'LOCALBV')
#' hpda_spe_cell_dec
#' gsvaExp(hpda_spe_cell_dec, 'LOCALBV')
#' # Then using ggsc to visualize the result
#' \donttest{
#'   library(ggplot2)
#'   library(ggsc)
#'   gsvaExp(hpda_spe_cell_dec, 'LOCALBV') |>
#'   plot_lisa_feature(res1, assay.type='LocalLee') + ggtitle(NULL)
#' }
setGeneric('runLOCALBV',
  function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    bv.method = c("locallee", "localmoran_bv"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    lisa.flag.method = c("mean", "median"),
    reduction.used = NULL,
    group.by = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  )
  standardGeneric('runLOCALBV')
)


#' @rdname runLOCALBV-method
#' @aliases runLOCALBV,SingleCellExperiment
#' @export runLOCALBV
setMethod("runLOCALBV", "SingleCellExperiment", function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    sample_id = 'all',
    bv.method = c("locallee", "localmoran"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    lisa.flag.method = c("mean", "median"),
    reduction.used = NULL,
    group.by = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  ){

  weight.method <- match.arg(weight.method)
  bv.method <- match.arg(bv.method)
  lisa.method <- match.arg(lisa.method)
  lisa.flag.method <- match.arg(lisa.flag.method)
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
  x <- .check_dgCMatrix(x)
  weight <- .check_weight(data, sample_id, weight, group.by)

  coords <- .check_coords(data, reduction.used, weight, weight.method)
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
                  wi <- Matrix::rowSums(wm)
                  wi2 <- Matrix::rowSums(wm^2)
                  if (any(Matrix::rowSums(wm) == 0)){
                      cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                  }
                  result <- .runLocalBv(xi, wm, features1, features2, n, NULL, across.gsvaexp, 
                                        permutation, bv.method, bv.alternative, random.seed, wi, 
                                        wi2, lisa.method, lisa.alternative, lisa.flag.method, BPPARAM)
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
    features1 = NULL,
    features2 = NULL,
    assay.type = 'logcounts',
    sample_id = 'all',
    bv.method = c("locallee", "localmoran_bv"),
    bv.alternative = "two.sided",
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    lisa.method = c("localG", "localmoran"),
    lisa.alternative = "greater",
    lisa.flag.method = c("mean", "median"),
    reduction.used = NULL,
    group.by = NULL,
    permutation = 100,
    random.seed = 1024,
    BPPARAM = SerialParam(),
    action = c("get", "only", "add"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  ){

    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform(c("The {.var gsvaexp} is specified, the specified {.var gsvaExp} will be",
                           "used to perform the analysis with `features` displayed from `rownames(data)`."))
       }
       weight.method <- match.arg(weight.method)
       bv.method <- match.arg(bv.method)
       lisa.method <- match.arg(lisa.method)
       lisa.flag.method <- match.arg(lisa.flag.method)
       action <- match.arg(action)
       sample_id <- .check_sample_id(data, sample_id)

       if (is.null(assay.type)){
           assay.type <- 1
       }

       x <- assay(data, assay.type)
       if (is.null(features1) && is.null(features2) && is.null(gsvaexp.features)){
           cli::cli_abort(c("The {.var gsvaexp} is specified, and the `data` is {.cls {class(data)}}.",
                            "The `features1`, `features2` and `gsvaexp.features` should not be `NULL` simultaneously."))
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
       x <- .check_dgCMatrix(x)
       listn <- .generate_feature_listn(data, features1, features2, gsvaexp)

       weight <- .check_weight(data, sample_id, weight, group.by)
       coords <- .check_coords(data, reduction.used, weight, weight.method)

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
                       wi <- Matrix::rowSums(wm)
                       wi2 <- Matrix::rowSums(wm^2)
                       if (any(Matrix::rowSums(wm) == 0)){
                           cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                       }
                       result <- .runLocalBv(xi, wm, features1, features2, n, listn, across.gsvaexp, permutation, bv.method, bv.alternative,
                                             random.seed, wi, wi2, lisa.method, lisa.alternative, lisa.flag.method, BPPARAM)
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
