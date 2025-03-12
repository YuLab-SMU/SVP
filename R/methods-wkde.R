#' @title Calculating the 2D Weighted Kernel Density Estimation
#' 
#' @rdname runWKDE-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param reduction.used character used as spatial coordinates to detect SVG, default is NULL,
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates, if this is provided
#' the coordinate of specified reduced result will be used.
#' @param grid.n integer number of grid points in the two directions to estimate 2D weighted kernel density, default is 100.
#' @param adjust numeric to adjust the \code{bandwidths}, default is 1.0.
#' @param bandwidths vector a two length numeric vector represents the bandwidths for x and y directions, default is normal 
#' reference bandwidth \code{\link[ks]{hpi}}, see also \code{\link[MASS]{bandwidth.nrd}}.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... additional parameters
#' @return a \linkS4class{SVPExperiment} or \linkS4class{SingleCellExperiment}
#' @export
#' @author Shuangbin Xu
#' @examples
#' library(SpatialExperiment)
#' data(hpda_spe_cell_dec)
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |> runWKDE(assay.type = 'affi.score')
#' # The result is saved in the assays (affi.score.density name) of SVPExperiment
#' # which can be extracted using assay and visualized using ggsc or
#' # other packages
#' assays(hpda_spe_cell_dec)
#' #\donttest{
#'     library(ggsc)
#'     
#'     f1 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone A",
#'                mapping=aes(x=x,y=y),
#'                slot = 'affi.score.density',
#'                pointsize=10
#'     ) +
#'     scale_bg_color_manual(values=c('black'))
#'     f1
#'     
#'     f2 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone B",
#'                mapping=aes(x=x,y=y),
#'                pointsize=10,
#'                slot = 'affi.score.density'
#'     ) +
#'     scale_bg_color_manual(values=c('black'))
#'     f2
#' #}
setGeneric('runWKDE',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = NULL,
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  )
  standardGeneric('runWKDE')
)

#' @importFrom SummarizedExperiment assay<-
#' @rdname runWKDE-method
#' @aliases runWKDE,SingleCellExperiment
#' @export runWKDE
setMethod('runWKDE', 'SingleCellExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = NULL,
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,    
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){

  if (is.null(assay.type)){
      assay.type <- assayNames(data)[1]
  }
 
  if (is.numeric(assay.type)){
      assay.type <- assayNames(data)[assay.type]
  } 

  x <- assay(data, assay.type)

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(reduction.used %in% reducedDimNames(data))

  if(flag1 || flag2){
      if (flag1){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }

      if (flag2){
          coords <- reducedDim(data, reduction.used)
          coords <- coords[,c(1, 2)]
      }      
      
      tic()
      if (verbose){
          cli::cli_inform("Running the 2D Weighted Density Estimation using RcppParallel ...")
      }
      new.assay.nm <- paste0(assay.type, ".density")
      res <- .internal_runWKDE(x, coords, bandwidths, adjust, grid.n)
      
      assay(data, new.assay.nm) <- res
      toc()
  }else{
      cli::cli_abort("The {.cls {class(data)}} should have 'spatialCoords' or the reduction result of 'UMAP' or 'TSNE'.")
  }
  return(data)
})

##' @rdname runWKDE-method   
##' @aliases runWKDE,SVPExperiment   
##' @export runWKDE
setMethod('runWKDE', 'SVPExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = NULL,
    grid.n = 100,
    adjust = 1,
    bandwidths = NULL,    
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){

    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to calculate density.")
       }
       da2 <- gsvaExp(data, gsvaexp, withSpatialCoords=TRUE, withReducedDim=TRUE, withColData = FALSE, withImgData = FALSE)
       da2 <- runWKDE(da2, 
                   assay.type = gsvaexp.assay.type, 
                   reduction.used = reduction.used, 
                   grid.n = grid.n, 
                   adjust = adjust, 
                   bandwidths = bandwidths)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
})

.internal_runWKDE <- function(x, coords, bandwidths = NULL, adjust = 1, grid.n = 100){

  lims <- c(range(coords[,1]), range(coords[,2]))

  if (is.null(bandwidths)){
      bandwidths <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
  }
  x <- .check_dgCMatrix(x)
  res <- CalWkdeParallel(x=as.matrix(coords), w=x, l=lims, h = bandwidths, adjust = adjust, n = grid.n)
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  return(res)
}
