#' @title Local indicators of spatial association analysis 
#' @description
#' This function use the local indicators of spatial association (LISA) to identify the hotspot
#' in the spatial space.
#' @rdname runLISA-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param features the feature name or index of data object, which are required.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param method character one of \code{'localG'}, \code{"localmoran"}, default is \code{'localG'}.
#' @param weight object, which can be \code{nb}, \code{listw} or \code{Graph} object, default is NULL.
#' @param weight.method character the method to build the spatial neighbours weights, default
#' is \code{knn} (k nearest neighbours). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is avaiable, such as \code{'tri2nb'}, \code{"knearneigh"},
#' \code{'dnearneigh'}, \code{"gabrielneigh"}, \code{"relativeneigh"}, which are from \code{spdep} package.
#' default is \code{knn}, if it is \code{"none"}, meaning the distance weight of each spot is used to
#' the weight.
#' @param reduction.used character used as spatial coordinates to calculate the neighbours weights, 
#' default is \code{UMAP}, if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param cells the cell name or index of data object, default is NULL.
#' @param action character, which control the type of return result, default is \code{get}, which will return
#' a \linkS4class{SimpleList}.
#' @param alternative a character string specifying the alternative hypothesis, default is \code{two.sided}.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis parallelly using
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel.
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)} to parallel it,
#' the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}..
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return if \code{action = 'get'} (in default), the SimpleList object (like list object) will be return, 
#' if \code{action = 'only'}, the data.frame will be return. if \code{action = 'add'}, the result of LISA is 
#' stored in the \code{localResults} column of \code{int_colData} (internal column metadata). You can use 
#' \code{localResults()} function of \code{SpatialFeatureExperiment} package to extract it.
#' @references
#' 1. Bivand, R.S., Wong, D.W.S. Comparing implementations of global and local indicators of spatial association. TEST 27,
#'    716–748 (2018). https://doi.org/10.1007/s11749-018-0599-x
#' @export
setGeneric('runLISA',
  function(
    data,
    features,
    assay.type = 'logcounts',
    method = c("localG", "localmoran"),
    weight = NULL,
    weight.method = c("knn", "tri2nb"),
    reduction.used = NULL,
    cells = NULL,
    action = c("get", "add", "only"),
    alternative = 'two.sided',
    BPPARAM = SerialParam(),
    ...
  )
  standardGeneric('runLISA')
)

#' @rdname runLISA-method
#' @aliases runLISA,SingleCellExperiment
#' @export runLISA
setMethod("runLISA", "SingleCellExperiment", function(
    data, 
    features, 
    assay.type = "logcounts",
    method = c("localG", "localmoran"),
    weight = NULL, 
    weight.method = c("knn", "tri2nb"), 
    reduction.used = NULL,
    cells = NULL,
    action = c("get", "add", "only"),
    alternative = 'two.sided',
    BPPARAM = SerialParam(),
    ...
  ){
  
  action <- match.arg(action)

  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)

  if (!is.null(cells)){
      x <- x[, cells, drop=FALSE]
  }

  x <- x[features, ,drop=FALSE]

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(reduction.used %in% reducedDimNames(data))


  if((flag1 || flag2)){
      if (flag2){
          coords <- reducedDim(data, reduction.used)
          coords <- coords[,c(1, 2)]
      }
      if (flag1 ){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
  }else if (all(flag1, flag2) && is.null(weight)){
      cli::cli_abort("The {.cls {class(data)}} should have 'spatialCoords' or the reduction result of 'UMAP' or 'TSNE'.
                     Or the `weight` should be provided.")
  }

  wm <- .obtain.weight(coords, weight = weight, weight.method = weight.method, ...)

  n <- nrow(wm)
  wi <- rowSums(wm)
  wi2 <- rowSums(wm^2)
  if (any(wi == 0)){
      cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
  }

  res <- BiocParallel::bplapply(seq(nrow(x)), function(i){ 
                                .internal.runLISA(x[i,], wm, wi, wi2, n, method, alternative)
             }, BPPARAM = BPPARAM           
         )

  names(res) <- rownames(x)

  if (action == 'only'){
      res <- do.call("cbind", res)
      return(res)
  }

  if (action == 'add'){
      data <- .add_LISA_res(data, res, method)
      return(data)
  }
  res <- S4Vectors::SimpleList(res)
  return(res)
})

