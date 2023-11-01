#' @title Calculating the 2D Weighted Kernel Density Estimation
#' @rdname wkde-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param reduction character used as spatial coordinates to detect SVG, default is \code{UMAP},
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param grid.n integer number of grid points in the two directions to estimate 2D weighted kernel density, default is 100.
#' @param adjust numeric to adjust the \code{bandwidths}, default is 1.0.
#' @param bandwidths vector a two length numeric vector represents the bandwidths for x and y directions, default is normal 
#' reference bandwidth \code{\link[ks]{hpi}}, see also \code{\link[MASS]{bandwidth.nrd}}.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... additional parameters
#' @export
setGeneric('wkde',
  function(
    data,
    assay.type = 'logcounts',
    reduction = c('UMAP', 'TSNE'),
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  )
  standardGeneric('wkde')
)

#' @importFrom SummarizedExperiment assay<-
#' @rdname wkde-method
#' @aliases wkde,SingleCellExperiment
#' @export wkde
setMethod('wkde', 'SingleCellExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction = c('UMAP', 'TSNE'),
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,    
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){
  if (is.numeric(assay.type)){
      assay.type <- assayNames(data)[assay.type]
  }      
  if (!assay.type %in% assayNames(data)){
      cli::cli_abort("the {.var assay.type} = {assay.type} is not present in the assays of {.cls {class(data)}}.")
  }

  x <- assay(data, assay.type)

  flag <- .check_element_obj(data, key = 'gsvaExps', basefun=int_colData, namefun = names)
  if (!is.null(gsvaexp) && flag){
    data <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withImgData = FALSE, withReducedDim = TRUE)
    if (is.numeric(gsvaexp.assay.type)){gsvaexp.assay.type <- assayNames(data)[gsvaexp.assay.type]}
    if (is.null(gsvaexp.assay.type)){gsvaexp.assay.type <- assayNames(data)[1]}
    x <- assay(data, gsvaexp.assay.type)
    assay.type <- gsvaexp.assay.type
  }

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(reduction %in% reducedDimNames(data))

  if(flag1 || flag2){
      if (flag2){
          coords <- reducedDim(data, reduction)
          coords <- coords[,c(1, 2)]
      }
      if (flag1){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
      
      tic()
      if (verbose){
          cli::cli_inform("Running the 2D Weighted Density Estimation using RcppParallel ...")
      }
      new.assay.nm <- paste0(assay.type, ".density")
      res <- .internal_wkde(x, coords, bandwidths, adjust, grid.n)
      
      assay(data, new.assay.nm) <- res
      toc()
  }else{
      cli::cli_abort("The {.cls {class(data)}} should have 'spatialCoords' or the reduction result of 'UMAP' or 'TSNE'.")
  }
  return(data)
})

##' @rdname wkde-method	  
##' @aliases wkde,SVPExperiment	  
##' @export wkde
setMethod('wkde', 'SVPExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction = c('UMAP', 'TSNE'),
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,    
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){

    x <- callNextMethod()
    #if (!is.null(gsvaexp)){
    #   if (verbose){
    #      cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to calculate density.")
    #   }
    #   gsvaExp(data, gsvaexp) <- x
    #   return(data)
    #}
    return(x)
})

.internal_wkde <- function(x, coords, bandwidths = NULL, adjust = 1, grid.n = 100){
  coords <- .normalize.coords(coords)

  lims <- c(range(coords[,1]), range(coords[,2]))

  if (is.null(bandwidths)){
      bandwidths <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
  }

  res <- CalWkdeParallel(x=as.matrix(coords), w=x, l=lims, h = bandwidths, adjust = adjust, n = grid.n)
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  return(res)
}	    
